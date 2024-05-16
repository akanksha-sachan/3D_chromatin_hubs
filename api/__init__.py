# __init__.py for web-package (API)

print("Initializing genomic-hubs package...")

__all__ = ["submodule1", "submodule2"] #define which namespaces/modules from the codebase are imported when the package is imported

# initialization code
def _init_package():
    print("Package genomic-hubs initialized.")

_init_package()
