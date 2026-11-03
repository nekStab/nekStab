"""Each LAPACK wrapper must return the result of rank 0 on all ranks.

Intel MKL can give different bits on different MPI ranks. A wrapper that returns
its local result lets the ranks rotate the Krylov basis differently and take
different decisions, which hangs the run in a global sum. See the header of
src/lapack_wrapper.f90.
"""

import re
from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parents[2] / "src" / "lapack_wrapper.f90"

OUTPUTS = {
    "schur": ["A", "vecs", "vals"],
    "ordschur": ["T", "Q"],
    "eig": ["vals", "vecs"],
    "lstsq": ["x"],
    "eig_symmetric": ["eigvals", "eigvecs"],
    "eig_hermitian": ["eigvals", "eigvecs"],
}


def body(name: str) -> str:
    text = SOURCE.read_text()
    match = re.search(rf"subroutine {name}\(.*?end subroutine {name}\b", text, re.S)
    assert match, f"{name} not found"
    return match.group(0)


@pytest.mark.parametrize("name, outputs", OUTPUTS.items())
def test_wrapper_broadcasts_its_outputs(name, outputs):
    code = body(name)
    for out in outputs:
        assert re.search(rf"call bcast\({out},", code), f"{name} does not broadcast {out}"


def test_every_wrapper_is_listed():
    text = SOURCE.read_text()
    lapack_calls = re.findall(r"call (dgees|dtrsen|dgeev|dgels|dsyev|zheev)\(", text)
    assert set(lapack_calls) == {"dgees", "dtrsen", "dgeev", "dgels", "dsyev", "zheev"}
    wrappers = set(re.findall(r"^\s*subroutine (\w+)\(", text, re.M)) - {"sort_eigendecomp"}
    assert wrappers == set(OUTPUTS)
