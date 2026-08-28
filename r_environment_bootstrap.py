import os
import sys
from pathlib import Path


def configure_r_environment() -> None:
    """
    Prefer the active conda environment's R installation when present.
    This must run before importing rpy2.
    """
    sys_prefix = Path(sys.prefix)
    env_r_home = sys_prefix / "lib" / "R"
    env_r_lib = env_r_home / "lib"
    env_r_library = env_r_home / "library"

    if not env_r_home.exists():
        return

    os.environ.setdefault("R_HOME", str(env_r_home))
    os.environ.setdefault("R_LIBS", str(env_r_library))
    os.environ.setdefault("R_LIBS_USER", str(env_r_library))

    ld_parts = [str(env_r_lib)]
    existing_ld = os.environ.get("LD_LIBRARY_PATH")
    if existing_ld:
        ld_parts.append(existing_ld)
    os.environ["LD_LIBRARY_PATH"] = ":".join(part for part in ld_parts if part)


configure_r_environment()
