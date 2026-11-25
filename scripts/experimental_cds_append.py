"""Compatibility entry; prefer `pantrans append`."""
from pantrans.cds_append import *  # Backward-compatible helper imports.
from pantrans.cds_append import main

if __name__ == "__main__":
    raise SystemExit(main())
