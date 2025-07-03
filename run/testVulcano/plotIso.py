import pyvista as pv
import glob
import os
import re
import time
import numpy as np
import argparse
import shutil
import subprocess
from multiprocessing import Pool, cpu_count

# ----- DEFAULT SETTINGS -----
DEFAULT_FRAMERATE = 10
DEFAULT_VIDEO_RESOLUTION = (1920, 1080)
DEFAULT_NUM_PROCESSES = max(1, cpu_count() - 1)
BACKGROUND_COLOR = 'white'

# ----- ARGUMENT PARSER -----
parser = argparse.ArgumentParser(
    description="""
Generate videos from OpenFOAM surface outputs.
... (description as before) ...
""",
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(
    "--interactive",
    action="store_true",
    help=
    "Enable interactive view to rotate texture, save camera (for isometric base), and then proceed."
)
parser.add_argument(
    "--np",
    "--num_processes",
    type=int,
    default=None,
    help=
    f"Number of parallel processes for frame generation (default: {DEFAULT_NUM_PROCESSES} or auto-detected)."
)
parser.add_argument(
    "--fr",
    "--framerate",
    type=int,
    default=None,
    help=f"Framerate for the output videos (default: {DEFAULT_FRAMERATE}).")
args = parser.parse_args()

# ----- SETTINGS FROM ARGS OR DEFAULTS -----
FRAMERATE = args.fr if args.fr is not None else DEFAULT_FRAMERATE
NUM_PROCESSES = args.np if args.np is not None else DEFAULT_NUM_PROCESSES
NUM_PROCESSES = max(1, min(NUM_PROCESSES, cpu_count()))
VIDEO_RESOLUTION = DEFAULT_VIDEO_RESOLUTION

# ----- GLOBAL VARIABLES FOR CAMERA -----
continue_animation_flag = False
user_defined_isometric_base_view = None
calculated_top_down_view = None

# ----- HELPER FUNCTIONS -----
# ... (read_foam_dict_value come prima) ...


def read_foam_dict_value(filepath, keyword):
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        content = re.sub(r'//.*', '', content)
        content = re.sub(r'/\*.*?\*/', '', content, flags=re.DOTALL)
        match = re.search(rf'{keyword}\s+([\d.eE+-]+);', content)
        if match:
            return float(match.group(1))
        else:
            print(f"Warning: Keyword '{keyword}' not found in '{filepath}'.")
            return None
    except FileNotFoundError:
        print(f"Warning: Dictionary file '{filepath}' not found.")
        return None
    except Exception as e:
        print(f"Error reading '{keyword}' from '{filepath}': {e}")
        return None


def key_press_callback_interactive(plotter_instance):
    global continue_animation_flag, user_defined_isometric_base_view
    user_defined_isometric_base_view = plotter_instance.camera_position
    print(
        f"User-defined isometric base view saved: {user_defined_isometric_base_view}"
    )
    continue_animation_flag = True
    plotter_instance.close()


def rotate_view_around_focal_point(camera_view, angle_deg):
    pos, focal, view_up = camera_view
    angle_rad = np.deg2rad(angle_deg)
    vec_focal_to_pos = np.array(pos) - np.array(focal)
    R_z = np.array([[np.cos(angle_rad), -np.sin(angle_rad), 0],
                    [np.sin(angle_rad),
                     np.cos(angle_rad), 0], [0, 0, 1]])
    vec_focal_to_pos_rotated = R_z @ vec_focal_to_pos
    new_pos = np.array(focal) + vec_focal_to_pos_rotated
    new_view_up = R_z @ np.array(view_up)
    return (new_pos.tolist(), focal, new_view_up.tolist())


# ----- PART 1: TERRAIN/TEXTURE LOADING & CAMERA SETUP -----
TERRAIN_VTK_PATH = "postProcessing/terrain/0/terrain.vtk"
TEXTURE_IMAGE_PATH = "texture.jpg"  # Stringa del percorso

print("Loading terrain mesh...")
if not os.path.exists(TERRAIN_VTK_PATH):
    raise FileNotFoundError(
        f"FATAL: Terrain file not found: {TERRAIN_VTK_PATH}")
terrain_mesh = pv.read(TERRAIN_VTK_PATH)
if not terrain_mesh or terrain_mesh.n_points == 0:
    raise ValueError(
        f"FATAL: Terrain mesh from '{TERRAIN_VTK_PATH}' is empty or invalid.")
print(
    f"Terrain mesh loaded. Bounds: {terrain_mesh.bounds}, Center: {terrain_mesh.center}"
)
terrain_mesh.compute_normals(auto_orient_normals=True,
                             consistent_normals=True,
                             inplace=True)
terrain_mesh.texture_map_to_plane(inplace=True)

texture_exists_globally = os.path.exists(TEXTURE_IMAGE_PATH)  # Flag globale
terrain_color_fallback = "lightgray"
# Non carichiamo texture_object globalmente se deve essere passato ai worker
# Lo caricheremo solo per la visualizzazione interattiva e i worker lo caricheranno da sé.

if texture_exists_globally:
    print(
        f"Texture file '{TEXTURE_IMAGE_PATH}' found. It will be loaded by workers or interactive plotter."
    )
else:
    print(
        f"Warning: Texture file '{TEXTURE_IMAGE_PATH}' not found. Using solid color."
    )

texture_rotation_angle_deg = [0]
texture_flip_x = [False]
texture_flip_y = [False]


def apply_texture_transformations_to_coords(mesh_with_uvs):
    if "Texture Coordinates" not in mesh_with_uvs.point_data:
        mesh_with_uvs.texture_map_to_plane(inplace=True)
    t_coords = mesh_with_uvs.active_t_coords.copy()
    t_coords[:, 0] -= 0.5
    t_coords[:, 1] -= 0.5
    angle_rad = np.deg2rad(texture_rotation_angle_deg[0])
    u, v = t_coords[:, 0].copy(), t_coords[:, 1].copy()
    t_coords[:, 0] = np.cos(angle_rad) * u - np.sin(angle_rad) * v
    t_coords[:, 1] = np.sin(angle_rad) * u + np.cos(angle_rad) * v
    t_coords[:, 0] += 0.5
    t_coords[:, 1] += 0.5
    if texture_flip_x[0]:
        t_coords[:, 0] = 1.0 - t_coords[:, 0]
    if texture_flip_y[0]:
        t_coords[:, 1] = 1.0 - t_coords[:, 1]
    mesh_with_uvs.active_t_coords = t_coords
    return mesh_with_uvs


if texture_exists_globally:  # Applica trasformazioni UV a terrain_mesh globalmente
    terrain_mesh = apply_texture_transformations_to_coords(terrain_mesh)

# --- Camera Setup Logic ---


def get_specific_view(mesh, view_type='iso'):
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(mesh.copy())
    if view_type == 'iso':
        plotter.camera_position = 'iso'
        plotter.reset_camera(render=False)
    elif view_type == 'xy':
        plotter.view_xy(negative=False)
    cam_pos = plotter.camera_position
    plotter.close()
    del plotter
    return cam_pos


interactive_texture_object = None  # Solo per la sessione interattiva
if args.interactive:
    print("\nStarting interactive setup for ISOMETRIC BASE VIEW...")
    if texture_exists_globally:  # Carica la texture solo per la finestra interattiva
        interactive_texture_object = pv.read_texture(TEXTURE_IMAGE_PATH)

    plotter_interactive = pv.Plotter(window_size=VIDEO_RESOLUTION,
                                     off_screen=False)
    plotter_interactive.background_color = BACKGROUND_COLOR

    actor_terrain_interactive = None
    if texture_exists_globally and interactive_texture_object:
        actor_terrain_interactive = plotter_interactive.add_mesh(
            terrain_mesh,
            texture=interactive_texture_object,
            smooth_shading=True)
    else:
        actor_terrain_interactive = plotter_interactive.add_mesh(
            terrain_mesh, color=terrain_color_fallback, smooth_shading=True)
    plotter_interactive.add_axes(zlabel='Z (Up)', xlabel='X', ylabel='Y')

    initial_iso_cam = get_specific_view(terrain_mesh, 'iso')
    plotter_interactive.camera_position = initial_iso_cam
    print(
        f"Initial camera for interactive mode (default iso): {plotter_interactive.camera_position}"
    )

    def cb_rotate_texture_interactive():
        texture_rotation_angle_deg[0] = (texture_rotation_angle_deg[0] +
                                         90) % 360
        print(f"Texture rotation: {texture_rotation_angle_deg[0]} degrees")
        global terrain_mesh
        # Start with a fresh copy of original geometry
        terrain_mesh_copy = terrain_mesh.copy()
        terrain_mesh_copy.texture_map_to_plane(
            inplace=True)  # Reset UVs to default
        terrain_mesh = apply_texture_transformations_to_coords(
            terrain_mesh_copy)  # Apply new transforms
        # L'actor usa il terrain_mesh globale che ora ha le TCoords aggiornate
        plotter_interactive.update_scalars(terrain_mesh.active_t_coords,
                                           mesh=actor_terrain_interactive,
                                           render=True)

    def cb_flip_texture_x_interactive():
        texture_flip_x[0] = not texture_flip_x[0]
        print(f"Texture flip X: {texture_flip_x[0]}")
        global terrain_mesh
        terrain_mesh_copy = terrain_mesh.copy()
        terrain_mesh_copy.texture_map_to_plane(inplace=True)
        terrain_mesh = apply_texture_transformations_to_coords(
            terrain_mesh_copy)
        plotter_interactive.update_scalars(terrain_mesh.active_t_coords,
                                           mesh=actor_terrain_interactive,
                                           render=True)

    def cb_flip_texture_y_interactive():
        texture_flip_y[0] = not texture_flip_y[0]
        print(f"Texture flip Y: {texture_flip_y[0]}")
        global terrain_mesh
        terrain_mesh_copy = terrain_mesh.copy()
        terrain_mesh_copy.texture_map_to_plane(inplace=True)
        terrain_mesh = apply_texture_transformations_to_coords(
            terrain_mesh_copy)
        plotter_interactive.update_scalars(terrain_mesh.active_t_coords,
                                           mesh=actor_terrain_interactive,
                                           render=True)

    plotter_interactive.add_key_event("r", cb_rotate_texture_interactive)
    plotter_interactive.add_key_event("x", cb_flip_texture_x_interactive)
    plotter_interactive.add_key_event("y", cb_flip_texture_y_interactive)
    plotter_interactive.add_key_event(
        "Return", lambda: key_press_callback_interactive(plotter_interactive))

    print(
        "  Interactive Window: Adjust ISOMETRIC base view. R:Rotate Txt, X/Y:Flip Txt, Enter:Save Cam, Q:Quit."
    )
    plotter_interactive.show(title="Interactive Isometric Base View Setup")

    if not continue_animation_flag:
        print(
            "Interactive setup exited without saving. Using default calculated isometric base view."
        )
        user_defined_isometric_base_view = get_specific_view(
            terrain_mesh, 'iso')

    if plotter_interactive.renderer:
        plotter_interactive.close()
    del plotter_interactive
    if interactive_texture_object:
        del interactive_texture_object  # Pulisci texture interattiva

else:
    print(
        "\nNon-interactive mode. Using default calculated isometric base view."
    )
    user_defined_isometric_base_view = get_specific_view(terrain_mesh, 'iso')

print("Calculating top-down view...")
calculated_top_down_view = get_specific_view(terrain_mesh, 'xy')

print(
    f"--- User-defined/Default ISOMETRIC BASE view for animations: {user_defined_isometric_base_view} ---"
)
print(
    f"--- Calculated TOP-DOWN view for animations: {calculated_top_down_view} ---"
)

# ----- PART 2: PREPARING ANIMATION DATA -----
print("\nLocating VTK files for animation...")
ISO_SURFACES_DIR_PATTERN = "postProcessing/isoSurfaces/*"
ISO4_VTK_NAME = "iso4.vtk"
ISO6_VTK_NAME = "iso6.vtk"


def get_timestep_from_dir_path(dir_path):
    try:
        return float(os.path.basename(dir_path))
    except ValueError:
        return -1


timestep_dirs = sorted(glob.glob(ISO_SURFACES_DIR_PATTERN),
                       key=get_timestep_from_dir_path)
valid_timestep_dirs = [
    d for d in timestep_dirs if get_timestep_from_dir_path(d) >= 0
]
iso4_vtk_files, iso6_vtk_files, simulation_timesteps = [], [], []
for ts_dir in valid_timestep_dirs:
    iso4_path = os.path.join(ts_dir, ISO4_VTK_NAME)
    iso6_path = os.path.join(ts_dir, ISO6_VTK_NAME)
    if os.path.exists(iso4_path) and os.path.exists(iso6_path):
        iso4_vtk_files.append(iso4_path)
        iso6_vtk_files.append(iso6_path)
        simulation_timesteps.append(get_timestep_from_dir_path(ts_dir))
    else:
        print(f"Warning: Missing VTK in {ts_dir}. Skipping.")
if not iso4_vtk_files:
    raise FileNotFoundError(
        f"No complete VTK sets found in {ISO_SURFACES_DIR_PATTERN}.")
print(
    f"Found {len(iso4_vtk_files)} timesteps with both iso4 and iso6 VTK files."
)

OUTPUT_FRAMES_ROOT_DIR = "postProcessing/frames_png"
OUTPUT_VIDEOS_ROOT_DIR = "postProcessing/videos_mp4"
os.makedirs(OUTPUT_FRAMES_ROOT_DIR, exist_ok=True)
os.makedirs(OUTPUT_VIDEOS_ROOT_DIR, exist_ok=True)

# ----- FUNCTION FOR RENDERING A SINGLE FRAME (WORKER) -----


def render_single_frame_worker(
        frame_index,
        num_total_frames,
        iso4_filepath,
        iso6_filepath,
        current_sim_time,
        camera_view_for_frame,
        output_png_filepath,
        passed_terrain_mesh,  # Questo è il terrain_mesh con le TCoords già trasformate
        passed_texture_exists_flag,  # Flag per sapere se la texture dovrebbe esistere
        passed_texture_image_path,  # Percorso del file texture
        passed_terrain_color):

    frame_plotter = pv.Plotter(off_screen=True, window_size=VIDEO_RESOLUTION)
    frame_plotter.background_color = BACKGROUND_COLOR

    worker_texture_object = None
    if passed_texture_exists_flag:
        try:
            worker_texture_object = pv.read_texture(passed_texture_image_path)
        except Exception as e:
            if frame_index == 0:
                print(
                    f"[WORKER {os.getpid()}] Error loading texture {passed_texture_image_path} in worker: {e}"
                )
            # Continua senza texture se il caricamento fallisce
            passed_texture_exists_flag = False  # Sovrascrivi il flag per questo worker

    if passed_texture_exists_flag and worker_texture_object and passed_terrain_mesh:
        # passed_terrain_mesh ha già le coordinate UV corrette
        frame_plotter.add_mesh(passed_terrain_mesh,
                               texture=worker_texture_object,
                               smooth_shading=True,
                               name="terrain")
    elif passed_terrain_mesh:
        frame_plotter.add_mesh(passed_terrain_mesh,
                               color=passed_terrain_color,
                               smooth_shading=True,
                               name="terrain")

    try:
        mesh4 = pv.read(iso4_filepath)
        mesh6 = pv.read(iso6_filepath)
        if mesh4.n_points > 0:
            frame_plotter.add_mesh(mesh4,
                                   color="tan",
                                   opacity=0.7,
                                   smooth_shading=True,
                                   name="iso4")
        if mesh6.n_points > 0:
            frame_plotter.add_mesh(mesh6,
                                   color="lightblue",
                                   opacity=0.7,
                                   smooth_shading=True,
                                   name="iso6")
    except Exception as e:
        if frame_index == 0:
            print(
                f"[WORKER {os.getpid()}] Error loading VTK for {output_png_filepath}: {e}"
            )

    frame_plotter.add_text(f"Time: {current_sim_time:.2f}s",
                           font_size=15,
                           position="upper_right",
                           color="black",
                           shadow=True)
    frame_plotter.add_axes(zlabel='Z', xlabel='X', ylabel='Y')

    frame_plotter.camera_position = camera_view_for_frame

    try:
        frame_plotter.screenshot(output_png_filepath,
                                 transparent_background=False)
    except Exception as e:
        print(
            f"[WORKER {os.getpid()}] EXCEPTION during screenshot for {output_png_filepath}: {e}"
        )

    if worker_texture_object:
        del worker_texture_object  # Pulisci texture del worker
    frame_plotter.close()
    del frame_plotter
    return output_png_filepath


# ----- FUNCTION TO ASSEMBLE VIDEO -----


def assemble_video_from_frames(frames_input_dir,
                               video_output_path,
                               framerate_val,
                               delete_frames_after=True):
    print(
        f"Assembling video: '{video_output_path}' from '{frames_input_dir}' at {framerate_val} FPS."
    )
    png_filename_pattern = os.path.join(frames_input_dir, "frame_%04d.png")
    ffmpeg_command = [
        'ffmpeg', '-r',
        str(framerate_val), '-i', png_filename_pattern, '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p', '-vf', "pad=ceil(iw/2)*2:ceil(ih/2)*2", '-an',
        '-y', video_output_path
    ]
    try:
        subprocess.run(ffmpeg_command,
                       check=True,
                       capture_output=True,
                       text=True)
        print(f"Video successfully saved: {video_output_path}")
        if delete_frames_after:
            shutil.rmtree(frames_input_dir)
            print(f"Directory '{frames_input_dir}' removed.")
    except subprocess.CalledProcessError as e:
        print(
            f"ERROR: FFMPEG failed for '{video_output_path}'.\n  Command: {' '.join(e.cmd)}\n  Stderr:\n{e.stderr}"
        )
        print(f"PNG frames preserved in '{frames_input_dir}'.")
    except FileNotFoundError:
        print("ERROR: FFMPEG not found. Ensure it's installed and in PATH.")
        print(f"PNG frames preserved in '{frames_input_dir}'.")


# ----- MAIN FUNCTION TO GENERATE SEQUENCE -----
def create_animation_sequence(base_camera_view,
                              rotation_angle_degrees=None,
                              use_circular_orbit=False,
                              sequence_base_name="video_seq"):
    current_sequence_frames_dir = os.path.join(OUTPUT_FRAMES_ROOT_DIR,
                                               sequence_base_name)
    os.makedirs(current_sequence_frames_dir, exist_ok=True)
    print(
        f"\nPreparing frames for sequence: '{sequence_base_name}' in '{current_sequence_frames_dir}'"
    )

    num_animation_frames = len(iso4_vtk_files)
    if num_animation_frames == 0:
        print("No frames to render. Skipping.")
        return

    worker_tasks_list = []
    for i in range(num_animation_frames):
        current_camera_for_frame = base_camera_view
        if use_circular_orbit:
            orbit_angle_deg = 360.0 * i / num_animation_frames
            current_camera_for_frame = rotate_view_around_focal_point(
                base_camera_view, orbit_angle_deg)
        elif rotation_angle_degrees is not None:
            current_camera_for_frame = rotate_view_around_focal_point(
                base_camera_view, rotation_angle_degrees)

        output_png_filepath = os.path.join(current_sequence_frames_dir,
                                           f"frame_{i:04d}.png")
        task_arg_tuple = (
            i,
            num_animation_frames,
            iso4_vtk_files[i],
            iso6_vtk_files[i],
            simulation_timesteps[i],
            current_camera_for_frame,
            output_png_filepath,
            terrain_mesh,  # terrain_mesh ha già le TCoords finali
            texture_exists_globally,  # Flag globale
            TEXTURE_IMAGE_PATH,  # Percorso del file texture
            terrain_color_fallback if not texture_exists_globally else None)
        worker_tasks_list.append(task_arg_tuple)

    print(
        f"Rendering {num_animation_frames} frames for '{sequence_base_name}' using {NUM_PROCESSES} processes..."
    )
    rendering_start_time = time.time()
    with Pool(processes=NUM_PROCESSES) as frame_pool:
        results = frame_pool.starmap(render_single_frame_worker,
                                     worker_tasks_list)
    rendering_end_time = time.time()
    print(
        f"Frames for '{sequence_base_name}' rendered in {rendering_end_time - rendering_start_time:.2f}s."
    )

    rendered_files_count = sum(1 for r in results if r and os.path.exists(r))
    if rendered_files_count != num_animation_frames:
        print(
            f"Warning: Expected {num_animation_frames} frames, got {rendered_files_count} for '{sequence_base_name}'."
        )

    final_video_output_path = os.path.join(OUTPUT_VIDEOS_ROOT_DIR,
                                           f"{sequence_base_name}.mp4")
    assemble_video_from_frames(current_sequence_frames_dir,
                               final_video_output_path,
                               FRAMERATE,
                               delete_frames_after=True)


# ----- MAIN EXECUTION SCRIPT -----
if __name__ == "__main__":
    print(
        f"Script started. Using {NUM_PROCESSES} processes for rendering at {FRAMERATE} FPS."
    )
    print(
        f"Video resolution: {VIDEO_RESOLUTION[0]}x{VIDEO_RESOLUTION[1]}, Background: {BACKGROUND_COLOR}"
    )

    if not user_defined_isometric_base_view or not calculated_top_down_view:
        raise RuntimeError(
            "FATAL: Base camera views were not properly initialized.")

    create_animation_sequence(base_camera_view=calculated_top_down_view,
                              sequence_base_name="video_view_top_down")

    static_isometric_angles = [0, 90, 180, 270]
    for angle_val in static_isometric_angles:
        create_animation_sequence(
            base_camera_view=user_defined_isometric_base_view,
            rotation_angle_degrees=angle_val,
            sequence_base_name=f"video_view_iso_{angle_val}deg")

    create_animation_sequence(
        base_camera_view=user_defined_isometric_base_view,
        use_circular_orbit=True,
        sequence_base_name="video_orbit_iso_360deg")

    print("\nAll requested video sequences processed.")
    print(f"Final videos should be in: '{OUTPUT_VIDEOS_ROOT_DIR}'")
