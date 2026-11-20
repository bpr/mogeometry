# ManifoldSurfaceMesh, its element handles and navigators, and MeshData,
# against geometry-central: tests/data/surface/reference.txt is written by
# tests/data/surface/reference.cpp (`pixi run -e ref ref-data`), which
# builds geometry-central's ManifoldSurfaceMesh from the same face lists
# and records its connectivity and element queries. The other tests port
# geometry-central's test/src/halfedge_mesh_test.cpp (those that apply to
# a manifold mesh without geometry).

from std.collections import Dict, Set
from std.testing import (
    TestSuite,
    assert_equal,
    assert_false,
    assert_raises,
    assert_true,
)

from mogeometry import (
    INVALID_IND,
    MESH_CHECKS,
    BoundaryLoop,
    BoundaryLoopData,
    Corner,
    DisjointSets,
    Edge,
    EdgeData,
    Face,
    FaceData,
    Halfedge,
    ManifoldSurfaceMesh,
    MarkedDisjointSets,
    Vertex,
    VertexData,
)

comptime REFERENCE = "tests/data/surface/reference.txt"


# ===----------------------------------------------------------------------=== #
# Reading the reference data
# ===----------------------------------------------------------------------=== #


struct MeshRecord(Copyable, Movable):
    """One mesh of the reference data: its faces and geometry-central's
    arrays (or, for an invalid input, whether geometry-central threw)."""

    var name: String
    var valid: Bool
    var threw: Bool
    var faces: List[List[Int]]
    var arrays: Dict[String, List[Int]]

    def __init__(out self, name: String, valid: Bool, threw: Bool):
        self.name = name
        self.valid = valid
        self.threw = threw
        self.faces = List[List[Int]]()
        self.arrays = Dict[String, List[Int]]()

    def array(self, name: String) raises -> List[Int]:
        return self.arrays[name].copy()


def _records() raises -> List[MeshRecord]:
    var out = List[MeshRecord]()
    with open(REFERENCE, "r") as f:
        var lines = f.read().split("\n")
        var i = 0
        while i < len(lines):
            var words = lines[i].split()
            i += 1
            if len(words) == 0:
                continue
            if words[0] != "mesh" and words[0] != "invalid":
                raise Error("unexpected line: " + String(lines[i - 1]))
            var valid = words[0] == "mesh"
            var threw = (not valid) and atol(words[2]) == 1
            var rec = MeshRecord(String(words[1]), valid, threw)
            var faces_words = lines[i].split()
            i += 1
            var n_faces = atol(faces_words[1])
            for _ in range(n_faces):
                var w = lines[i].split()
                i += 1
                var face = List[Int]()
                for k in range(1, len(w)):
                    face.append(atol(w[k]))
                rec.faces.append(face^)
            while True:
                var w = lines[i].split()
                i += 1
                if w[0] == "end":
                    break
                var xs = List[Int]()
                for k in range(2, len(w)):
                    xs.append(atol(w[k]))
                if len(xs) != atol(w[1]):
                    raise Error("bad array length: " + String(w[0]))
                rec.arrays[String(w[0])] = xs^
            out.append(rec^)
    return out^


def _meshes() raises -> List[MeshRecord]:
    var out = List[MeshRecord]()
    for r in _records():
        if r.valid:
            out.append(r.copy())
    if len(out) == 0:
        raise Error("no meshes in " + REFERENCE)
    return out^


def _mesh(name: String) raises -> ManifoldSurfaceMesh:
    for r in _meshes():
        if r.name == name:
            return ManifoldSurfaceMesh(r.faces)
    raise Error("no reference mesh named " + name)


def _expect(rec: MeshRecord, array: String, i: Int, got: Int, want: Int) raises:
    if got != want:
        raise Error(
            String(
                rec.name,
                ": ",
                array,
                "[",
                i,
                "] is ",
                got,
                ", geometry-central's is ",
                want,
            )
        )


# ===----------------------------------------------------------------------=== #
# Against geometry-central
# ===----------------------------------------------------------------------=== #


def test_connectivity_matches_geometry_central() raises:
    for rec in _meshes():
        var mesh = ManifoldSurfaceMesh(rec.faces)
        var counts = rec.array("counts")
        assert_equal(mesh.n_vertices(), counts[0], rec.name)
        assert_equal(mesh.n_edges(), counts[1], rec.name)
        assert_equal(mesh.n_faces(), counts[2], rec.name)
        assert_equal(mesh.n_halfedges(), counts[3], rec.name)
        assert_equal(mesh.n_interior_halfedges(), counts[4], rec.name)
        assert_equal(mesh.n_boundary_loops(), counts[5], rec.name)

        var he_next = rec.array("he_next")
        var he_vertex = rec.array("he_vertex")
        var he_face = rec.array("he_face")
        var he_edge = rec.array("he_edge")
        var he_interior = rec.array("he_interior")
        for i in range(mesh.n_halfedges()):
            var h = mesh.halfedge(i)
            _expect(rec, "he_next", i, mesh.next(h).index(), he_next[i])
            _expect(rec, "he_vertex", i, mesh.vertex(h).index(), he_vertex[i])
            _expect(rec, "he_face", i, mesh.face(h).index(), he_face[i])
            _expect(rec, "he_edge", i, mesh.edge(h).index(), he_edge[i])
            _expect(
                rec,
                "he_interior",
                i,
                Int(mesh.is_interior(h)),
                he_interior[i],
            )

        var v_halfedge = rec.array("v_halfedge")
        var v_degree = rec.array("v_degree")
        var v_face_degree = rec.array("v_face_degree")
        var v_boundary = rec.array("v_boundary")
        for i in range(mesh.n_vertices()):
            var v = mesh.vertex(i)
            _expect(
                rec, "v_halfedge", i, mesh.halfedge(v).index(), v_halfedge[i]
            )
            _expect(rec, "v_degree", i, mesh.degree(v), v_degree[i])
            _expect(
                rec, "v_face_degree", i, mesh.face_degree(v), v_face_degree[i]
            )
            _expect(
                rec, "v_boundary", i, Int(mesh.is_boundary(v)), v_boundary[i]
            )

        var f_halfedge = rec.array("f_halfedge")
        var f_degree = rec.array("f_degree")
        for i in range(mesh.n_faces()):
            var f = mesh.face(i)
            _expect(
                rec, "f_halfedge", i, mesh.halfedge(f).index(), f_halfedge[i]
            )
            _expect(rec, "f_degree", i, mesh.degree(f), f_degree[i])

        var e_halfedge = rec.array("e_halfedge")
        var e_boundary = rec.array("e_boundary")
        for i in range(mesh.n_edges()):
            var e = mesh.edge(i)
            _expect(
                rec, "e_halfedge", i, mesh.halfedge(e).index(), e_halfedge[i]
            )
            _expect(
                rec, "e_boundary", i, Int(mesh.is_boundary(e)), e_boundary[i]
            )

        var bl_halfedge = rec.array("bl_halfedge")
        var bl_face = rec.array("bl_face")
        var bl_degree = rec.array("bl_degree")
        for i in range(mesh.n_boundary_loops()):
            var b = mesh.boundary_loop(i)
            _expect(
                rec, "bl_halfedge", i, mesh.halfedge(b).index(), bl_halfedge[i]
            )
            _expect(rec, "bl_face", i, mesh.as_face(b).index(), bl_face[i])
            _expect(rec, "bl_degree", i, mesh.degree(b), bl_degree[i])

        var props = rec.array("properties")
        assert_equal(mesh.euler_characteristic(), props[0], rec.name)
        assert_equal(mesh.genus(), props[1], rec.name)
        assert_equal(mesh.n_connected_components(), props[2], rec.name)
        assert_equal(Int(mesh.is_triangular()), props[3], rec.name)
        assert_equal(Int(mesh.has_boundary()), props[4], rec.name)
        assert_equal(mesh.n_interior_vertices(), props[5], rec.name)


def test_validate_connectivity() raises:
    # ValidateClosedMeshTest, ValidateBoundaryMeshTest
    for rec in _meshes():
        ManifoldSurfaceMesh(rec.faces).validate_connectivity()


def test_invalid_input_raises() raises:
    var n = 0
    for rec in _records():
        if rec.valid:
            continue
        assert_true(rec.threw, "geometry-central accepted " + rec.name)
        with assert_raises():
            _ = ManifoldSurfaceMesh(rec.faces)
        n += 1
    assert_equal(n, 8)


def test_degenerate_input_raises() raises:
    var empty = List[List[Int]]()
    with assert_raises(contains="unreferenced vertex 0"):
        _ = ManifoldSurfaceMesh(empty)
    var negative: List[List[Int]] = [[0, -1, 2]]
    with assert_raises(contains="negative vertex index"):
        _ = ManifoldSurfaceMesh(negative)


# ===----------------------------------------------------------------------=== #
# geometry-central's halfedge_mesh_test.cpp
# ===----------------------------------------------------------------------=== #


def test_iterate_elements() raises:
    # IterateVerticesTest ... IterateBoundaryLoopsTest
    for rec in _meshes():
        var mesh = ManifoldSurfaceMesh(rec.faces)
        var vs = Set[Vertex]()
        for v in mesh.vertices():
            vs.add(v)
        assert_equal(len(vs), mesh.n_vertices(), rec.name)
        var hs = Set[Halfedge]()
        for h in mesh.halfedges():
            hs.add(h)
        assert_equal(len(hs), mesh.n_halfedges(), rec.name)
        var interior = Set[Halfedge]()
        for h in mesh.interior_halfedges():
            assert_true(mesh.is_interior(h))
            interior.add(h)
        assert_equal(len(interior), mesh.n_interior_halfedges(), rec.name)
        var exterior = Set[Halfedge]()
        for h in mesh.exterior_halfedges():
            assert_false(mesh.is_interior(h))
            exterior.add(h)
        assert_equal(len(exterior), mesh.n_exterior_halfedges(), rec.name)
        var cs = Set[Corner]()
        for c in mesh.corners():
            cs.add(c)
        assert_equal(len(cs), mesh.n_corners(), rec.name)
        var es = Set[Edge]()
        for e in mesh.edges():
            es.add(e)
        assert_equal(len(es), mesh.n_edges(), rec.name)
        var fs = Set[Face]()
        for f in mesh.faces():
            fs.add(f)
        assert_equal(len(fs), mesh.n_faces(), rec.name)
        var bs = Set[BoundaryLoop]()
        for b in mesh.boundary_loops():
            bs.add(b)
        assert_equal(len(bs), mesh.n_boundary_loops(), rec.name)


def test_index_count() raises:
    # IndexCountTest: every index is below the index size, and the largest
    # reaches it (corners may fall one short)
    for rec in _meshes():
        var mesh = ManifoldSurfaceMesh(rec.faces)
        var top = -1
        for v in mesh.vertices():
            assert_true(v.index() < mesh.vertex_index_size())
            top = max(top, v.index())
        assert_equal(top + 1, mesh.vertex_index_size(), rec.name)
        top = -1
        for h in mesh.halfedges():
            assert_true(h.index() < mesh.halfedge_index_size())
            top = max(top, h.index())
        assert_equal(top + 1, mesh.halfedge_index_size(), rec.name)
        top = -1
        for c in mesh.corners():
            assert_true(c.index() < mesh.corner_index_size())
            top = max(top, c.index())
        assert_true(top + 2 >= mesh.corner_index_size(), rec.name)
        top = -1
        for e in mesh.edges():
            assert_true(e.index() < mesh.edge_index_size())
            top = max(top, e.index())
        assert_equal(top + 1, mesh.edge_index_size(), rec.name)
        top = -1
        for f in mesh.faces():
            assert_true(f.index() < mesh.face_index_size())
            top = max(top, f.index())
        assert_equal(top + 1, mesh.face_index_size(), rec.name)


def test_whole_mesh_queries() raises:
    # HasBoundaryTest, IsTriangularTest, EulerCharacteristicTest, GenusTest,
    # ConnectedComponentsTest, on the meshes of the test that are here
    var tet = _mesh("tet")
    assert_false(tet.has_boundary())
    assert_true(tet.is_triangular())
    assert_equal(tet.euler_characteristic(), 2)
    assert_equal(tet.genus(), 0)
    assert_equal(tet.n_connected_components(), 1)
    var dodecahedron = _mesh("dodecahedron")
    assert_false(dodecahedron.is_triangular())
    assert_equal(dodecahedron.euler_characteristic(), 2)
    assert_equal(dodecahedron.genus(), 0)
    var bob = _mesh("bob_small")
    assert_true(bob.is_triangular())
    assert_equal(bob.euler_characteristic(), 0)
    assert_equal(bob.genus(), 1)
    assert_equal(bob.n_connected_components(), 1)
    assert_true(_mesh("grid_triangles").has_boundary())
    assert_equal(_mesh("shelf").n_connected_components(), 5)
    assert_true(tet.is_manifold())
    assert_true(tet.is_edge_manifold())
    assert_true(tet.is_oriented())


def test_print_element() raises:
    # PrintElementTest, PrintElementStringTest
    var mesh = _mesh("grid_quads")
    assert_equal(String(mesh.halfedge(0)), "he_0")
    assert_equal(String(mesh.corner(0)), "c_0")
    assert_equal(String(mesh.vertex(0)), "v_0")
    assert_equal(String(mesh.edge(3)), "e_3")
    assert_equal(String(mesh.face(0)), "f_0")
    assert_equal(String(mesh.boundary_loop(0)), "bl_0")
    assert_equal(String(Vertex()), "v_-1")


def test_prev_navigators() raises:
    # PrevTest
    for rec in _meshes():
        var mesh = ManifoldSurfaceMesh(rec.faces)
        for h in mesh.halfedges():
            var n = mesh.next(h)
            assert_true(mesh.prev_orbit_face(n) == h, rec.name)
            assert_true(mesh.prev_orbit_vertex(n) == h, rec.name)


def test_vertex_navigators() raises:
    # VertexAdjacentNavigator, VertexCornerNavigatorInterior,
    # VertexEdgeNavigator, VertexFaceNavigator; also incoming halfedges and
    # neighbors
    for rec in _meshes():
        var mesh = ManifoldSurfaceMesh(rec.faces)
        var outgoing_total = 0
        var edge_count = EdgeData[Int](mesh, 0)
        var face_count = FaceData[Int](mesh, 0)
        for v in mesh.vertices():
            for h in mesh.outgoing_halfedges(v):
                assert_true(mesh.vertex(h) == v)
                outgoing_total += 1
            var n_in = 0
            for h in mesh.incoming_halfedges(v):
                assert_true(mesh.tip_vertex(h) == v)
                n_in += 1
            assert_equal(n_in, mesh.degree(v), rec.name)
            var n_adj = 0
            for w in mesh.adjacent_vertices(v):
                assert_true(mesh.connecting_edge(v, w).is_valid())
                n_adj += 1
            assert_equal(n_adj, mesh.degree(v), rec.name)
            for c in mesh.adjacent_corners(v):
                assert_false(mesh.is_boundary_loop(mesh.face(c)))
                assert_true(mesh.vertex(c) == v)
            var seen_edges = Set[Edge]()
            for e in mesh.adjacent_edges(v):
                assert_false(e in seen_edges, rec.name)
                seen_edges.add(e)
                edge_count[e] += 1
            var seen_faces = Set[Face]()
            for f in mesh.adjacent_faces(v):
                assert_false(f in seen_faces, rec.name)
                seen_faces.add(f)
                face_count[f] += 1
        assert_equal(outgoing_total, mesh.n_halfedges(), rec.name)
        for e in mesh.edges():
            assert_equal(edge_count[e], 2, rec.name)
        for f in mesh.faces():
            assert_equal(face_count[f], mesh.degree(f), rec.name)


def test_face_navigators() raises:
    # FaceFaceNavigator; also the halfedge, vertex, corner, and edge walks
    for rec in _meshes():
        var mesh = ManifoldSurfaceMesh(rec.faces)
        var face_count = FaceData[Int](mesh, 0)
        for f in mesh.faces():
            var seen = Set[Face]()
            for g in mesh.adjacent_faces(f):
                assert_false(g in seen, rec.name)
                seen.add(g)
                face_count[g] += 1
            var k = 0
            var vertices = List[Vertex]()
            for v in mesh.adjacent_vertices(f):
                vertices.append(v)
            for h in mesh.adjacent_halfedges(f):
                assert_true(mesh.face(h) == f)
                assert_true(mesh.vertex(h) == vertices[k])
                k += 1
            assert_equal(k, mesh.degree(f))
            for c in mesh.adjacent_corners(f):
                assert_true(mesh.face(c) == f)
            for e in mesh.adjacent_edges(f):
                var touches = False
                for g in mesh.adjacent_faces(e):
                    touches = touches or g == f
                assert_true(touches)
        for f in mesh.faces():
            var expected = 0
            for e in mesh.adjacent_edges(f):
                expected += mesh.degree(e) - 1
            assert_equal(face_count[f], expected, rec.name)


def test_face_vertices_are_the_input() raises:
    # getFaceVertexList gives back the faces the mesh was built from
    for rec in _meshes():
        var mesh = ManifoldSurfaceMesh(rec.faces)
        var faces = mesh.get_face_vertex_list()
        assert_equal(len(faces), len(rec.faces))
        for i in range(len(faces)):
            assert_true(faces[i] == rec.faces[i], rec.name)


def test_edge_navigators() raises:
    # EdgeVertexNavigator, EdgeDiamondNavigator
    for rec in _meshes():
        var mesh = ManifoldSurfaceMesh(rec.faces)
        for e in mesh.edges():
            var ends = mesh.adjacent_vertices(e)
            assert_true(ends.first != ends.second)
            assert_true(ends.first == mesh.first_vertex(e))
            assert_true(ends.second == mesh.second_vertex(e))
            assert_true(mesh.other_vertex(e, ends.first) == ends.second)
            assert_true(mesh.other_vertex(e, ends.second) == ends.first)
            assert_true(mesh.connecting_edge(ends.first, ends.second) == e)
            var n = 0
            for h in mesh.adjacent_halfedges(e):
                assert_true(mesh.edge(h) == e)
                n += 1
            assert_equal(n, 2)
            assert_equal(mesh.degree(e), 1 if mesh.is_boundary(e) else 2)
    for name in ["tet", "bob_small", "torus", "bipyramid"]:
        var mesh = _mesh(name)
        for e in mesh.edges():
            var h = mesh.halfedge(e)
            var used = Set[Halfedge]()
            used.add(h)
            used.add(mesh.twin(h))
            var d = mesh.diamond_boundary(e)
            for k in range(4):
                assert_false(d[k] in used)
                used.add(d[k])
                assert_true(
                    mesh.face(d[k]) == mesh.face(h)
                    or mesh.face(d[k]) == mesh.face(mesh.twin(h))
                )
    var grid = _mesh("grid_triangles")
    for e in grid.edges():
        if grid.is_boundary(e):
            with assert_raises(contains="boundary edge"):
                _ = grid.diamond_boundary(e)
            break
    var quads = _mesh("grid_quads")
    for e in quads.edges():
        if not quads.is_boundary(e):
            with assert_raises(contains="non-triangular"):
                _ = quads.diamond_boundary(e)
            break
    var far = _mesh("tet")
    assert_false(far.connecting_edge(far.vertex(0), far.vertex(0)).is_valid())
    var fan = _mesh("fan")
    assert_false(fan.connecting_edge(fan.vertex(1), fan.vertex(10)).is_valid())


def test_boundary_loops() raises:
    for rec in _meshes():
        var mesh = ManifoldSurfaceMesh(rec.faces)
        var exterior = 0
        for b in mesh.boundary_loops():
            var f = mesh.as_face(b)
            assert_true(mesh.is_boundary_loop(f))
            assert_true(mesh.as_boundary_loop(f) == b)
            var k = 0
            for h in mesh.adjacent_halfedges(b):
                assert_false(mesh.is_interior(h))
                assert_true(mesh.face(h) == f)
                k += 1
            var vs = 0
            for v in mesh.adjacent_vertices(b):
                assert_true(mesh.is_boundary(v))
                vs += 1
            var es = 0
            for e in mesh.adjacent_edges(b):
                assert_true(mesh.is_boundary(e))
                es += 1
            assert_equal(k, mesh.degree(b))
            assert_equal(vs, k)
            assert_equal(es, k)
            exterior += k
        assert_equal(exterior, mesh.n_exterior_halfedges(), rec.name)
        for f in mesh.faces():
            assert_false(mesh.is_boundary_loop(f))
        for v in mesh.vertices():
            if mesh.is_boundary(v):
                # the interior halfedge along the boundary
                assert_true(mesh.is_interior(mesh.halfedge(v)))
                assert_false(mesh.is_interior(mesh.twin(mesh.halfedge(v))))
                assert_equal(mesh.face_degree(v), mesh.degree(v) - 1)
            else:
                assert_equal(mesh.face_degree(v), mesh.degree(v))


def test_halfedge_queries() raises:
    var mesh = _mesh("perforated")
    for h in mesh.halfedges():
        var t = mesh.twin(h)
        assert_true(mesh.twin(t) == h)
        assert_true(mesh.sibling(h) == t)
        assert_true(mesh.edge(t) == mesh.edge(h))
        assert_true(mesh.tip_vertex(h) == mesh.tail_vertex(t))
        assert_true(mesh.orientation(h) != mesh.orientation(t))
        assert_true(mesh.orientation(h) == (mesh.halfedge(mesh.edge(h)) == h))
        assert_true(mesh.next_outgoing_neighbor(h) == mesh.next(mesh.twin(h)))
        assert_true(mesh.next_incoming_neighbor(h) == mesh.twin(mesh.next(h)))
        assert_false(mesh.is_dead(h))
        assert_true(mesh.halfedge(mesh.corner(h)) == h)
    for v in mesh.vertices():
        assert_false(mesh.is_dead(v))
        assert_true(mesh.vertex(mesh.corner(v)) == v)
        assert_true(mesh.is_manifold(v))
    for e in mesh.edges():
        assert_false(mesh.is_dead(e))
        assert_true(mesh.is_manifold(e) and mesh.is_oriented(e))
    for f in mesh.faces():
        assert_false(mesh.is_dead(f))
        assert_false(mesh.is_triangle(f))
    for b in mesh.boundary_loops():
        assert_false(mesh.is_dead(b))
    assert_true(mesh.is_compressed())


def test_index_queries() raises:
    var mesh = _mesh("annulus")
    var vi = mesh.get_vertex_indices()
    for v in mesh.vertices():
        assert_equal(vi[v], v.index())
    var hi = mesh.get_halfedge_indices()
    for h in mesh.halfedges():
        assert_equal(hi[h], h.index())
    var ci = mesh.get_corner_indices()
    var k = 0
    for c in mesh.corners():
        assert_equal(ci[c], k)
        k += 1
    var ei = mesh.get_edge_indices()
    for e in mesh.edges():
        assert_equal(ei[e], e.index())
    var fi = mesh.get_face_indices()
    for f in mesh.faces():
        assert_equal(fi[f], f.index())
    var bi = mesh.get_boundary_loop_indices()
    for b in mesh.boundary_loops():
        assert_equal(bi[b], b.index())
    var interior = mesh.get_interior_vertex_indices()
    var boundary = mesh.get_vertex_boundary_status()
    var n = 0
    for v in mesh.vertices():
        if boundary[v]:
            assert_equal(interior[v], INVALID_IND)
        else:
            assert_equal(interior[v], n)
            n += 1
    assert_equal(n, mesh.n_interior_vertices())


# ===----------------------------------------------------------------------=== #
# Containers
# ===----------------------------------------------------------------------=== #


def test_container_access() raises:
    # ContainerAccessTest (with a loop over raw() for `raw() += raw()`)
    var mesh = _mesh("bob_small")
    var a = VertexData[Float64](mesh, 0.0)
    var b = VertexData[Float64](mesh, 0.0)
    for v in mesh.vertices():
        a[v] = 3.0
        b[v] = 2.0
    for i in range(len(a)):
        a.raw()[i] += b.raw()[i]
    for v in mesh.vertices():
        assert_equal(a[v], 5.0)
    assert_equal(len(a), mesh.n_vertices())
    assert_equal(a[4], 5.0)
    a.fill(1.5)
    assert_equal(a[mesh.vertex(7)], 1.5)
    assert_equal(a.default(), 0.0)
    a.set_default(9.0)
    assert_equal(a.default(), 9.0)


def test_container_sizes() raises:
    var mesh = _mesh("perforated")
    assert_equal(len(VertexData[Int](mesh, 0)), mesh.n_vertices())
    assert_equal(len(EdgeData[Int](mesh, 0)), mesh.n_edges())
    assert_equal(len(FaceData[Int](mesh, 0)), mesh.n_faces())
    assert_equal(len(BoundaryLoopData[Int](mesh, 0)), mesh.n_boundary_loops())
    var values = List[Float64]()
    for i in range(mesh.n_faces()):
        values.append(Float64(i))
    var data = FaceData[Float64](mesh, values^, 0.0)
    for f in mesh.faces():
        assert_equal(data[f], Float64(f.index()))
    with assert_raises(contains="values for"):
        _ = FaceData[Float64](mesh, List[Float64](), 0.0)


def test_container_sync() raises:
    var mesh = _mesh("torus")
    var data = VertexData[Int](mesh, 7)
    data.sync(mesh)
    assert_equal(len(data), mesh.n_vertices())
    var copy = mesh.copy()
    data.sync(copy)
    var other = _mesh("annulus")
    with assert_raises(contains="not the mesh"):
        data.sync(other)


# ===----------------------------------------------------------------------=== #
# Handles and mesh tags
# ===----------------------------------------------------------------------=== #


def test_handles() raises:
    var mesh = _mesh("tet")
    var v = mesh.vertex(2)
    assert_true(v == Vertex(2))
    assert_true(v.is_valid())
    assert_false(Vertex().is_valid())
    assert_equal(Vertex().index(), INVALID_IND)
    assert_true(Vertex(1) < Vertex(2))
    assert_true(Vertex(3) >= v)
    assert_equal(hash(v), hash(Vertex(2)))
    if MESH_CHECKS:
        assert_equal(v.tag(), mesh.mesh_tag())
        assert_equal(Vertex(2).tag(), 0)
    else:
        assert_equal(v.tag(), 0)


def test_owns() raises:
    var tet = _mesh("tet")
    var torus = _mesh("torus")
    assert_true(tet.owns(tet.vertex(3)))
    assert_true(tet.owns(Vertex(3)))  # untagged
    assert_false(tet.owns(Vertex(4)))  # out of range
    assert_false(tet.owns(Vertex()))
    var far = torus.vertex(3)
    if MESH_CHECKS:
        assert_true(tet.mesh_tag() != torus.mesh_tag())
        assert_false(tet.owns(far))
        assert_false(tet.owns(torus.halfedge(0)))
    else:
        assert_true(tet.owns(far))
    # a copy has the same elements, and the same tag
    var copy = tet.copy()
    assert_equal(copy.mesh_tag(), tet.mesh_tag())
    assert_true(copy.owns(tet.vertex(3)))
    assert_true(copy.next(tet.halfedge(5)) == tet.next(tet.halfedge(5)))


def test_mesh_ids_differ() raises:
    var tags = Set[Int]()
    for rec in _meshes():
        tags.add(ManifoldSurfaceMesh(rec.faces).mesh_tag())
    assert_equal(len(tags), len(_meshes()))


# ===----------------------------------------------------------------------=== #
# Disjoint sets
# ===----------------------------------------------------------------------=== #


def test_disjoint_sets() raises:
    var sets = DisjointSets(9)
    sets.merge(0, 1)
    sets.merge(2, 3)
    sets.merge(1, 3)
    sets.merge(5, 9)
    assert_equal(sets.find(0), sets.find(2))
    assert_equal(sets.find(5), sets.find(9))
    assert_true(sets.find(4) != sets.find(0))
    assert_true(sets.find(5) != sets.find(0))
    assert_equal(sets.find(4), 4)
    var marked = MarkedDisjointSets(5)
    marked.mark(1)
    assert_true(marked.is_marked(1))
    assert_false(marked.is_marked(2))
    marked.merge(2, 1)
    assert_true(marked.is_marked(2))
    marked.unmark(1)
    assert_false(marked.is_marked(2))
    assert_false(marked.is_marked(4))


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
