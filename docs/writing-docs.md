# Writing doc strings

The API reference will be generated from the doc strings (as MojoDEC's
is), so they are the documentation. Every public module, struct, trait, field,
alias, and function needs one.

## Format

Google style, as Mojo's standard library uses. The first paragraph is the
summary; then any longer description, in Markdown; then the sections
that apply:

```mojo
def compute_distance(mut self, source: Vertex) raises -> VertexData[Float64]:
    """The geodesic distance from `source` to every vertex, by the heat
    method.

    Distances are shifted so that the source's is 0.

    Args:
        source: The vertex the distance is measured from.

    Returns:
        The distance at each vertex of the mesh.

    Raises:
        If `source` is not a vertex of this mesh.
    """
```

- **`Parameters:`** documents compile-time parameters, on the struct or
  function that declares them (`MeshData[E, T]`).
- **`Args:`** lists every argument except `self`.
- **`Returns:`** is needed whenever the function returns a value.
- **`Raises:`** is needed whenever the function is `raises`.
- Fields and `comptime` aliases take a one-line doc string below them.
- A module's doc string starts with its title, underlined with `=`; the
  generator uses the title in the index.

Code in doc strings is not compiled or run (Mojo has no doctest runner),
so keep examples short and check them by hand; the tests and the examples
are the checked usage.

## Checking

```sh
pixi run doc-check                       # fail if any file got worse
pixi run doc-check --list surface/heat   # the warnings in matching files
pixi run doc-check --update              # record the counts after fixing some
```

`mojo doc --diagnose-missing-doc-strings` warns about missing doc strings
and missing sections. `doc-check` counts its warnings per file and fails
when a file has more than `docs/docstring-baseline.txt` records. Every
file is at zero, so the baseline is empty and any warning fails the
check: document new code as you write it. (`--update` exists for a
future tool version that adds new kinds of warning; it should not be
used to excuse undocumented code.)

## Ported code

A function ported from geometry-central says so in its doc string only
when the behavior differs (a different default, an error where
geometry-central returns garbage); the file header names the
geometry-central file it follows (see AGENTS.md, Porting).
