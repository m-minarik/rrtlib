# -*- coding: utf-8 -*-
"""
Creates an animation of an object following a single path from a specififed run.

Written for Blender 3.4
2023, Michal Minarik
"""
import sys
import bpy  # type: ignore
import json
import math

# useful shortcuts
scene = bpy.context.scene
collection = bpy.context.collection


def import_off(filepath: str, color=None):
    """
    Parse mesh from OFF file

    Supports COFF - colored off files
    x y z r g b a (rgba in [0, 1])

    If the input file is a OFF file without colors, a single RGBA color can be
    specified by the color parameter, which will be used for all vertices
    """
    file = open(filepath, "r")
    first_line = file.readline().rstrip()
    use_colors = first_line == "COFF"
    colors = []

    # handle blank and comment lines after the first line
    line = file.readline()
    while line.isspace() or line[0] == "#":
        line = file.readline()

    vcount, fcount, ecount = [int(x) for x in line.split()]
    verts = []
    facets = []
    edges = []
    i = 0
    while i < vcount:
        line = file.readline()
        if line.isspace():
            continue  # skip empty lines
        try:
            bits = [float(x) for x in line.split()]
            px = bits[0]
            py = bits[1]
            pz = bits[2]
            if use_colors:
                # r g b a in [0, 1]
                colors.append(bits[3:7])

        except ValueError:
            i = i + 1
            continue
        verts.append((px, py, pz))
        i = i + 1

    i = 0
    while i < fcount:
        line = file.readline()
        if line.isspace():
            continue  # skip empty lines
        try:
            splitted = line.split()
            ids = list(map(int, splitted))
            if len(ids) > 3:
                facets.append(tuple(ids[1:]))
            elif len(ids) == 3:
                edges.append(tuple(ids[1:]))
        except ValueError:
            i = i + 1
            continue
        i = i + 1

    # Assemble mesh
    off_name = bpy.path.display_name_from_filepath(filepath)
    mesh = bpy.data.meshes.new(name=off_name)
    mesh.from_pydata(verts, edges, facets)

    mesh.validate()
    mesh.update()

    color_data = mesh.vertex_colors.new()
    for i, facet in enumerate(mesh.polygons):
        for j, vidx in enumerate(facet.vertices):
            if use_colors:
                color_data.data[3 * i + j].color = colors[vidx]
            elif color is not None:
                color_data.data[3 * i + j].color = color
            else:
                color_data.data[3 * i + j].color = (0.9, 0.9, 0.9, 1.0)

    # Add new material
    mat = bpy.data.materials.new(name="Material")
    mesh.materials.append(mat)

    # Enable 'Use nodes' and add Vertex Color Node
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    vc = nodes.new("ShaderNodeVertexColor")
    vc.layer_name = "Col"

    # Get the shader
    bsdf = mat.node_tree.nodes["Principled BSDF"]

    # Link the vertex color and alpha to the shader
    mat.node_tree.links.new(vc.outputs[0], bsdf.inputs[0])
    mat.node_tree.links.new(vc.outputs[1], bsdf.inputs[18])

    # return mesh
    object = bpy.data.objects.new(off_name, mesh)
    collection.objects.link(object)
    return object


def cylinder_between(x1, y1, z1, x2, y2, z2, r, mat=None):
    """
    https://blender.stackexchange.com/questions/5898/how-can-i-create-a-cylinder-linking-two-points-with-python
    """
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    dist = math.sqrt(dx**2 + dy**2 + dz**2)

    if dist == 0:
        return

    bpy.ops.mesh.primitive_cylinder_add(
        vertices=8, radius=r, depth=dist, location=(dx / 2 + x1, dy / 2 + y1, dz / 2 + z1)
    )

    phi = math.atan2(dy, dx)
    theta = math.acos(dz / dist)

    bpy.context.object.rotation_euler[1] = theta
    bpy.context.object.rotation_euler[2] = phi

    if mat is not None:
        bpy.context.active_object.data.materials.append(mat)


def create_animation(run_dir, path_id=0):
    # Clear everything
    scene.camera = None
    for obj in bpy.context.visible_objects:
        bpy.context.collection.objects.unlink(obj)

    # Set light
    lamp_data = bpy.data.lights.new(name="lamp", type="POINT")
    lamp_object = bpy.data.objects.new(name="Lamp", object_data=lamp_data)
    bpy.context.collection.objects.link(lamp_object)
    lamp_object.location = (4, 1, 6)
    lamp = bpy.data.lights[lamp_data.name]
    lamp.energy = 1000
    lamp.use_shadow = False
    lamp.shadow_soft_size = 0.1

    # Set camera
    cam_data = bpy.data.cameras.new(name="camera")
    cam_ob = bpy.data.objects.new(name="Camera", object_data=cam_data)
    bpy.context.collection.objects.link(cam_ob)
    scene.camera = cam_ob
    cam_ob.location = (20, 10, 14.5)
    cam_ob.rotation_euler = (1.012291, 0, 2.042035)
    cam = bpy.data.cameras[cam_data.name]
    cam.lens = 50

    # Create material
    mat_name = "Trajectory"
    mat = bpy.data.materials.new(name=mat_name)
    mat.use_nodes = True
    mat.node_tree.nodes["Principled BSDF"].inputs[0].default_value = (0, 0.01, 0.8, 1)

    # Load object and map data
    with open(run_dir + "info.json", mode="r") as f:
        info = json.load(f)

    _ = import_off(info["params"]["map"])
    object = import_off(info["params"]["object"], color=(0, 0.01, 0.8, 1))

    # Create animation
    frame_no = 0
    p1 = None
    p2 = None
    s = 5

    with open(info["paths"][f"{path_id}"]["path_filename"], mode="r") as f:
        # Number of steps
        total_steps = int(f.readline().strip())

        for line in f:
            scene.frame_set(frame_no)

            # Set position and rotation (using quaternions)
            try:
                x, y, z, qx, qy, qz, qw = list(map(float, line.strip().split()))
            except Exception as e:
                continue

            object.rotation_mode = "QUATERNION"
            object.rotation_quaternion = (qw, qx, qy, qz)
            object.keyframe_insert(data_path="rotation_quaternion", index=-1)

            # Position
            object.location = (x, y, z)
            object.keyframe_insert(data_path="location", index=-1)

            frame_no += 1

            if s < 5:
                s += 1
                continue
            s = 0

            p1 = p2
            p2 = list(map(float, line.strip().split()))[0:3]
            if p1 is not None:
                cylinder_between(*p1, *p2, 0.02, mat)

    scene.frame_start = 1
    scene.frame_end = frame_no
    scene.frame_step = 1


create_animation(sys.argv[-1])
