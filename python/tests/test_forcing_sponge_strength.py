from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[2]


def _forcing_source() -> str:
    return (REPO_ROOT / "src" / "forcing.f90").read_text()


def _source(path: str) -> str:
    return (REPO_ROOT / path).read_text()


def _compact(source: str) -> str:
    return re.sub(r"\s+", "", source.lower())


def test_all_sponge_forcing_branches_scale_by_strength():
    source = _compact(_forcing_source())

    expected_terms = [
        "ffx=ffx+spng_fn(ip)*(spng_vr(ip,1)-vx(ix,iy,iz,iel))*spng_st",
        "ffy=ffy+spng_fn(ip)*(spng_vr(ip,2)-vy(ix,iy,iz,iel))*spng_st",
        "ffz=ffz+spng_fn(ip)*(spng_vr(ip,ndim)-vz(ix,iy,iz,iel))*spng_st",
        "ffx=ffx-spng_fn(ip)*vxp(ip,jp)*spng_st",
        "ffy=ffy-spng_fn(ip)*vyp(ip,jp)*spng_st",
        "ffz=ffz-spng_fn(ip)*vzp(ip,jp)*spng_st",
        "temp=temp+spng_fn(ip)*(spng_vt(ip,m)-t(ix,iy,iz,iel,m))*spng_st",
        "temp=temp-spng_fn(ip)*tp(ip,m,jp)*spng_st",
    ]

    for term in expected_terms:
        assert term in source


def test_sponge_mass_array_covers_temperature_elements():
    main = _compact(_source("src/main.f90"))
    copy_pos = main.find("callcopy(bm1s,bm1,nt)")

    assert "bm1s(lx1, ly1, lz1, lelt)" in _source("src/NEKSTAB")
    assert copy_pos != -1
    assert "nt=nx1*ny1*nz1*nelt" in main[:copy_pos]
