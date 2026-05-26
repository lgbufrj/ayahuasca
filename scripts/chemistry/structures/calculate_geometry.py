"""Compute transfer geometry (distance and angle) between donor and acceptor
atoms given pdb-mapped entries produced by the pipeline.

Functions
- distance(a, b): Euclidean distance
- angle_degrees(a, b, c): angle at b (in degrees)
- compute_transfer_geometry(pdb_mapped, transfer_vector): returns dict with
  donor/acceptor coords, distance, and angle (if a vertex is provided)

A small CLI at the bottom allows quick testing using a pipeline JSON file.
"""
from typing import Dict, Any, Optional
import math


def _distance(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _angle_degrees(a, b, c):
    # angle at b
    ba = (a[0] - b[0], a[1] - b[1], a[2] - b[2])
    bc = (c[0] - b[0], c[1] - b[1], c[2] - b[2])
    dot = ba[0] * bc[0] + ba[1] * bc[1] + ba[2] * bc[2]
    na = math.sqrt(ba[0] ** 2 + ba[1] ** 2 + ba[2] ** 2)
    nb = math.sqrt(bc[0] ** 2 + bc[1] ** 2 + bc[2] ** 2)
    if na == 0 or nb == 0:
        return None
    cosv = max(-1.0, min(1.0, dot / (na * nb)))
    return math.degrees(math.acos(cosv))


def compute_transfer_geometry(pdb_mapped: Dict[str, Any], transfer_vector: Dict[str, Any]) -> Dict[str, Any]:
    """Compute donor-acceptor geometry.

    Parameters
    - pdb_mapped: dict mapping atom_map (int or str) -> mapping entries produced
      by pdb_coords.align_atommap_and_compute_distances (entry contains 'coord').
    - transfer_vector: dict with 'donor' and 'acceptor' keys referencing atom_map ints.

    Returns
    - dict with donor, acceptor (each has map, coord, element), distance, and
      optionally angle if a vertex atom_map is provided in transfer_vector.
    """
    out = {
        'donor': None,
        'acceptor': None,
        'distance': None,
        'angle_at_acceptor': None,
    }

    try:
        dmap = transfer_vector.get('donor', {}).get('atom_map')
        amap = transfer_vector.get('acceptor', {}).get('atom_map')
    except Exception:
        return out

    if dmap is None or amap is None:
        return out

    # keys in pdb_mapped may be strings in JSON outputs
    def get_entry(k):
        return pdb_mapped.get(k) or pdb_mapped.get(str(k))

    dentry = get_entry(dmap)
    aentry = get_entry(amap)

    if not dentry or not aentry:
        return out

    dcoord = dentry.get('coord')
    acoord = aentry.get('coord')

    out['donor'] = {'map': dmap, 'coord': dcoord, 'element': dentry.get('element')}
    out['acceptor'] = {'map': amap, 'coord': acoord, 'element': aentry.get('element')}

    if dcoord and acoord:
        out['distance'] = round(_distance(dcoord, acoord), 3)

    # If a vertex is provided (third atom), compute angle at that atom
    vertex_map = None
    if 'vertex' in transfer_vector:
        vertex_map = transfer_vector.get('vertex')
    elif transfer_vector.get('acceptor', {}).get('vertex'):
        vertex_map = transfer_vector.get('acceptor').get('vertex')

    if vertex_map is not None:
        ventry = get_entry(vertex_map)
        if ventry and ventry.get('coord') and acoord:
            # compute angle donor-acceptor-vertex (or donor-vertex-acceptor?)
            # choose donor-vertex-acceptor (angle at vertex)
            vcoord = ventry.get('coord')
            # angle at vertex v: donor - vertex - acceptor
            ang = _angle_degrees(dcoord, vcoord, acoord)
            out['angle_at_vertex'] = round(ang, 2) if ang is not None else None

    return out


if __name__ == '__main__':
    # small smoke test: load a JSON output and compute transfer geometry
    import json
    import sys
    if len(sys.argv) < 2:
        print('Usage: python calculate_geometry.py /path/to/pipeline_output.json')
        sys.exit(1)
    jpath = sys.argv[1]
    data = json.load(open(jpath))
    pdb_mapped = data.get('pdb_mapping', {}).get('pdb_mapped', {})
    transfer_vector = data.get('transfer_vector', {})
    res = compute_transfer_geometry(pdb_mapped, transfer_vector)
    print(json.dumps(res, indent=2))
