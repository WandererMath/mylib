import os

def get_paths_ends_with_something(folder_path, suffix):
    return [os.path.abspath(os.path.join(folder_path, f)) for f in os.listdir(folder_path) if f.endswith(suffix)]
