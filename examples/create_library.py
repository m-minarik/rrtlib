import os
import pathlib

# Get absolute path to the root 'paper-rrtlib-code' directory
project_root = pathlib.Path(__file__).parent.parent.resolve()

# Set data directory paths
map_dir = f"{project_root}/data/maps"
obj_dir = f"{project_root}/data/objects"

# Prepare output directories and clear existing data
library = f"{project_root}/data/library/"
os.makedirs(library, exist_ok=True)
os.system(f"rm -rf {library}/*")

#     object      guiding_object    map     kwargs
runs = [
    ["chair112", "chair112_small", "1w3h", "-sc '3.5 3.0 2.0 1 0 0 0' -gc '-3.5 -3.0 2.0 1 0 0 0'"],
    ["table141", "table141_small", "1w3h", "-sc '3.5 3.0 2.0 1 0 0 0' -gc '-3.5 -3.0 2.0 1 0 0 0'"],
    ["teddy161", "teddy161_small", "1w3h", "-sc '3.5 3.0 2.0 1 0 0 0' -gc '-3.5 -3.0 2.0 1 0 0 0'"],
]

for r in runs:
    object, guiding_object, map, kwargs = r
    out_dir = f"{library}/{map}/{object}"
    os.makedirs(out_dir)

    # Prepare script for blender visualization
    blender_show = f"blender {project_root}/examples/blender/empty.blend -P {project_root}/examples/blender/visualize_trajectories.py -- {out_dir}/"
    with open(f"{out_dir}/blender_show.sh", mode="w") as f:
        f.write(blender_show)

    # Generate paths and save them to library
    job = f'{project_root}/main -M 2 -s 1 -m "{map_dir}/{map}.off" -o "{obj_dir}/{object}.off" \
        -go "{obj_dir}/{guiding_object}.off" --out "{out_dir}" {kwargs}'

    print(job)
