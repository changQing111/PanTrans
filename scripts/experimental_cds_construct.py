"""Compatibility entry; prefer `pantrans construct`."""
from pantrans.cds_construct import *  # Backward-compatible helper imports.
from pantrans.cds_construct import main

if __name__ == "__main__":
    raise SystemExit(main())
