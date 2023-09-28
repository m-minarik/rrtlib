import os
import pathlib

# Get absolute path to the root 'paper-rrtlib-code' directory
project_root = pathlib.Path(__file__).parent.parent.resolve()

# Set data directory paths
map_dir = f"{project_root}/data/maps"
obj_dir = f"{project_root}/data/objects"
library = f"{project_root}/data/library/"

object = f"{obj_dir}/chair111.off"
map = f"{map_dir}/1w3h.off"
kwargs = "-sc '3.5 3.0 2.0 1 0 0 0' -gc '-3.5 -3.0 2.0 1 0 0 0'"

out_dir = f"{project_root}/output/example"
os.makedirs(out_dir, exist_ok=True)
os.system(f"rm -rf {out_dir}/*")

# Identify the most similar object in the object library and save correspondences to output directory
corr_dir = f"{out_dir}/correspondences"
os.makedirs(corr_dir, exist_ok=True)
os.system(
    f"{project_root}/external/shape-similarity/identify-object -o '{object}' -l '{library}/1w3h' -t '{corr_dir}/temp' --out '{corr_dir}'"
)

# indentify-object will output correspondences and create symlinks to guiding_object and guiding_paths inside corr_dir
job = f"""
    {project_root}/main -M 1 -m '{map}' -o '{object}' \
    -go '{corr_dir}/guiding_object.off' -pd '{corr_dir}/guiding_paths/' \
    -c '{corr_dir}/correspondences.txt' --out '{out_dir}' \
    -s 1 {kwargs}
"""
os.system(job)

# Prepare script for blender visualization
blender_show = f"blender {project_root}/examples/blender/empty.blend -P {project_root}/examples/blender/create_animation.py -- {out_dir}/"
with open(f"{out_dir}/blender_show.sh", mode="w") as f:
    f.write(blender_show)