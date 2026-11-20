// Reference connectivity for tests/test_surface_mesh.mojo: geometry-central's
// ManifoldSurfaceMesh built from face lists, with its connectivity arrays
// and element queries. Build and run with `pixi run -e ref ref-data`; the
// first argument is the geometry-central checkout (for its CC0 test meshes).
//
// Output, one record per mesh:
//
//   mesh NAME
//   faces F
//   DEGREE V0 V1 ...           (F lines)
//   ARRAY N X0 X1 ...          (one line per array; see `arrays` below)
//   end
//
// and one per input geometry-central must reject:
//
//   invalid NAME THREW
//   faces F ... (as above)
//   end
//
// Every number is a decimal integer; -1 is geometry-central's INVALID_IND.

#include "geometrycentral/surface/manifold_surface_mesh.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

using Polygons = std::vector<std::vector<size_t>>;

static long long ind(size_t i) { return i == INVALID_IND ? -1 : (long long)i; }

static void putFaces(const Polygons& faces) {
  std::printf("faces %zu\n", faces.size());
  for (const auto& f : faces) {
    std::printf("%zu", f.size());
    for (size_t v : f) std::printf(" %zu", v);
    std::printf("\n");
  }
}

static void putArray(const char* name, const std::vector<long long>& xs) {
  std::printf("%s %zu", name, xs.size());
  for (long long x : xs) std::printf(" %lld", x);
  std::printf("\n");
}

static void emitMesh(const std::string& name, const Polygons& faces) {
  ManifoldSurfaceMesh mesh(faces);
  mesh.validateConnectivity();
  std::printf("mesh %s\n", name.c_str());
  putFaces(faces);

  std::vector<long long> heNext, heVertex, heFace, heEdge, heInterior;
  for (size_t i = 0; i < mesh.nHalfedges(); i++) {
    Halfedge he = mesh.halfedge(i);
    heNext.push_back(ind(he.next().getIndex()));
    heVertex.push_back(ind(he.vertex().getIndex()));
    heFace.push_back(ind(he.face().getIndex()));
    heEdge.push_back(ind(he.edge().getIndex()));
    heInterior.push_back(he.isInterior());
  }
  std::vector<long long> vHalfedge, vDegree, vFaceDegree, vBoundary;
  for (size_t i = 0; i < mesh.nVertices(); i++) {
    Vertex v = mesh.vertex(i);
    vHalfedge.push_back(ind(v.halfedge().getIndex()));
    vDegree.push_back(v.degree());
    vFaceDegree.push_back(v.faceDegree());
    vBoundary.push_back(v.isBoundary());
  }
  std::vector<long long> fHalfedge, fDegree;
  for (size_t i = 0; i < mesh.nFaces(); i++) {
    Face f = mesh.face(i);
    fHalfedge.push_back(ind(f.halfedge().getIndex()));
    fDegree.push_back(f.degree());
  }
  std::vector<long long> eHalfedge, eBoundary;
  for (size_t i = 0; i < mesh.nEdges(); i++) {
    Edge e = mesh.edge(i);
    eHalfedge.push_back(ind(e.halfedge().getIndex()));
    eBoundary.push_back(e.isBoundary());
  }
  std::vector<long long> blHalfedge, blFace, blDegree;
  for (size_t i = 0; i < mesh.nBoundaryLoops(); i++) {
    BoundaryLoop bl = mesh.boundaryLoop(i);
    blHalfedge.push_back(ind(bl.halfedge().getIndex()));
    blFace.push_back(ind(bl.asFace().getIndex()));
    blDegree.push_back(bl.degree());
  }

  putArray("counts", {(long long)mesh.nVertices(), (long long)mesh.nEdges(), (long long)mesh.nFaces(),
                      (long long)mesh.nHalfedges(), (long long)mesh.nInteriorHalfedges(),
                      (long long)mesh.nBoundaryLoops()});
  putArray("properties", {mesh.eulerCharacteristic(), mesh.genus(), (long long)mesh.nConnectedComponents(),
                          mesh.isTriangular(), mesh.hasBoundary(), (long long)mesh.nInteriorVertices()});
  putArray("he_next", heNext);
  putArray("he_vertex", heVertex);
  putArray("he_face", heFace);
  putArray("he_edge", heEdge);
  putArray("he_interior", heInterior);
  putArray("v_halfedge", vHalfedge);
  putArray("v_degree", vDegree);
  putArray("v_face_degree", vFaceDegree);
  putArray("v_boundary", vBoundary);
  putArray("f_halfedge", fHalfedge);
  putArray("f_degree", fDegree);
  putArray("e_halfedge", eHalfedge);
  putArray("e_boundary", eBoundary);
  putArray("bl_halfedge", blHalfedge);
  putArray("bl_face", blFace);
  putArray("bl_degree", blDegree);
  std::printf("end\n");
}

static void emitInvalid(const std::string& name, const Polygons& faces) {
  bool threw = false;
  try {
    ManifoldSurfaceMesh mesh(faces);
    mesh.validateConnectivity();
  } catch (const std::exception& e) {
    threw = true;
  }
  std::printf("invalid %s %d\n", name.c_str(), threw ? 1 : 0);
  putFaces(faces);
  std::printf("end\n");
}

// === Meshes from geometry-central's test assets (Keenan Crane's models, CC0)

// The face lists of a binary little-endian PLY with `property list uchar int`
// faces (the format of geometry-central's bob_small.ply and spot.ply)
static Polygons readPlyFaces(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::string line;
  size_t nVertices = 0, nFaces = 0, vertexBytes = 0;
  std::string element;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    std::string word;
    ss >> word;
    if (word == "element") {
      ss >> element;
      if (element == "vertex") ss >> nVertices;
      if (element == "face") ss >> nFaces;
    } else if (word == "property" && element == "vertex") {
      std::string type;
      ss >> type;
      if (type == "float" || type == "int" || type == "uint") vertexBytes += 4;
      else if (type == "double") vertexBytes += 8;
      else if (type == "uchar" || type == "char") vertexBytes += 1;
      else throw std::runtime_error("unexpected vertex property type " + type);
    } else if (word == "end_header") {
      break;
    }
  }
  in.seekg(nVertices * vertexBytes, std::ios::cur);
  Polygons faces(nFaces);
  for (size_t i = 0; i < nFaces; i++) {
    unsigned char degree;
    in.read(reinterpret_cast<char*>(&degree), 1);
    for (unsigned k = 0; k < degree; k++) {
      int32_t v;
      in.read(reinterpret_cast<char*>(&v), 4);
      faces[i].push_back(v);
    }
  }
  if (!in) throw std::runtime_error("truncated " + path);
  return faces;
}

// === Generated meshes

// tet.obj, 0-based
static Polygons tet() { return {{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2}}; }

// dodecahedron_poly.obj, 0-based: twelve pentagons
static Polygons dodecahedron() {
  Polygons f = {{4, 1, 2, 3, 5},     {8, 6, 7, 3, 2},     {12, 9, 10, 11, 13}, {16, 14, 15, 11, 10},
                {18, 17, 13, 11, 15}, {17, 18, 5, 3, 7},   {20, 19, 8, 2, 1},   {19, 20, 16, 10, 9},
                {19, 9, 12, 6, 8},    {6, 12, 13, 17, 7},  {14, 4, 5, 18, 15},  {4, 14, 16, 20, 1}};
  for (auto& p : f)
    for (auto& v : p) v -= 1;
  return f;
}

// An n x m grid of vertices; quads split into triangles along alternating
// diagonals if `triangles`, else kept as quads. A disk with one boundary loop.
static Polygons grid(size_t n, size_t m, bool triangles) {
  Polygons f;
  for (size_t i = 0; i + 1 < n; i++) {
    for (size_t j = 0; j + 1 < m; j++) {
      size_t a = i * m + j, b = a + 1, c = a + m + 1, d = a + m;
      if (!triangles) {
        f.push_back({a, b, c, d});
      } else if ((i + j) % 2 == 0) {
        f.push_back({a, b, c});
        f.push_back({a, c, d});
      } else {
        f.push_back({a, b, d});
        f.push_back({b, c, d});
      }
    }
  }
  return f;
}

// A quad grid with some quads left out, each hole away from the others and
// from the outer boundary: one boundary loop per hole plus the outer one.
static Polygons perforated(size_t n, std::vector<std::pair<size_t, size_t>> holes) {
  Polygons f;
  for (size_t i = 0; i + 1 < n; i++) {
    for (size_t j = 0; j + 1 < n; j++) {
      bool hole = false;
      for (auto h : holes) hole = hole || (h.first == i && h.second == j);
      if (hole) continue;
      size_t a = i * n + j;
      f.push_back({a, a + 1, a + n + 1, a + n});
    }
  }
  return f;
}

// k segments around, r rings of quads: two boundary loops
static Polygons annulus(size_t k, size_t r) {
  Polygons f;
  for (size_t i = 0; i < r; i++) {
    for (size_t j = 0; j < k; j++) {
      size_t a = i * k + j, b = i * k + (j + 1) % k;
      f.push_back({a, b, b + k, a + k});
    }
  }
  return f;
}

// An n x m torus, triangulated: closed, genus 1
static Polygons torus(size_t n, size_t m) {
  Polygons f;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      size_t a = i * m + j, b = i * m + (j + 1) % m;
      size_t c = ((i + 1) % n) * m + (j + 1) % m, d = ((i + 1) % n) * m + j;
      f.push_back({a, b, c});
      f.push_back({a, c, d});
    }
  }
  return f;
}

// Two cones of k triangles joined at their rims: closed, two vertices of
// degree k
static Polygons bipyramid(size_t k) {
  Polygons f;
  size_t top = k, bottom = k + 1;
  for (size_t j = 0; j < k; j++) {
    size_t a = j, b = (j + 1) % k;
    f.push_back({a, b, top});
    f.push_back({b, a, bottom});
  }
  return f;
}

// A fan of k triangles around vertex 0: a disk whose center has degree k
static Polygons fan(size_t k) {
  Polygons f;
  for (size_t j = 1; j <= k; j++) f.push_back({0, j, j % k + 1});
  return f;
}

// Disjoint union, renumbering each part's vertices after the previous ones
static Polygons combine(std::vector<Polygons> parts) {
  Polygons f;
  size_t offset = 0;
  for (auto& p : parts) {
    size_t nV = 0;
    for (auto& face : p) {
      std::vector<size_t> g;
      for (size_t v : face) {
        g.push_back(v + offset);
        nV = std::max(nV, v + 1);
      }
      f.push_back(g);
    }
    offset += nV;
  }
  return f;
}

// The same mesh with its vertices relabeled and its faces reordered and
// rotated at random (from a fixed seed)
static Polygons shuffled(Polygons f, uint64_t seed) {
  std::mt19937_64 rng(seed);
  size_t nV = 0;
  for (auto& p : f)
    for (size_t v : p) nV = std::max(nV, v + 1);
  std::vector<size_t> relabel(nV);
  for (size_t i = 0; i < nV; i++) relabel[i] = i;
  for (size_t i = nV - 1; i > 0; i--) std::swap(relabel[i], relabel[rng() % (i + 1)]);
  for (size_t i = f.size() - 1; i > 0; i--) std::swap(f[i], f[rng() % (i + 1)]);
  for (auto& p : f) {
    std::rotate(p.begin(), p.begin() + rng() % p.size(), p.end());
    for (auto& v : p) v = relabel[v];
  }
  return f;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s GEOMETRY_CENTRAL_PATH\n", argv[0]);
    return 1;
  }
  std::string assets = std::string(argv[1]) + "/test/assets/";

  emitMesh("tet", tet());
  emitMesh("dodecahedron", dodecahedron());
  emitMesh("bob_small", readPlyFaces(assets + "bob_small.ply"));
  emitMesh("grid_triangles", grid(6, 5, true));
  emitMesh("grid_quads", grid(5, 4, false));
  emitMesh("perforated", perforated(8, {{1, 1}, {1, 4}, {4, 2}, {5, 5}}));
  emitMesh("annulus", annulus(9, 3));
  emitMesh("torus", torus(8, 6));
  emitMesh("torus_shuffled", shuffled(torus(8, 6), 7));
  emitMesh("perforated_shuffled", shuffled(perforated(8, {{1, 1}, {1, 4}, {4, 2}, {5, 5}}), 11));
  emitMesh("bipyramid", bipyramid(40));
  emitMesh("fan", fan(30));
  emitMesh("shelf", combine({tet(), dodecahedron(), torus(4, 5), grid(4, 4, true), annulus(5, 1)}));

  // Inputs geometry-central rejects
  emitInvalid("hourglass", combine({tet(), {{0, 4, 5}, {0, 5, 6}, {0, 6, 4}, {4, 6, 5}}}));
  emitInvalid("bowtie", {{0, 1, 2}, {0, 3, 4}});
  emitInvalid("edge_three_faces", {{0, 1, 2}, {1, 0, 3}, {0, 1, 4}});
  emitInvalid("misoriented", {{0, 1, 2}, {0, 1, 3}});
  emitInvalid("duplicate_face", {{0, 1, 2}, {0, 1, 2}});
  emitInvalid("self_edge", {{0, 1, 1}, {0, 2, 1}});
  emitInvalid("degree_two", {{0, 1, 2}, {0, 1}});
  emitInvalid("unreferenced", {{0, 1, 2}, {0, 2, 4}});
  return 0;
}
