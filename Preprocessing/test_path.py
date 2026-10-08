import pathlib
import sys

def find_repo_root(start_path: pathlib.Path, repo_name: str) -> pathlib.Path:
    """
    Finds the absolute path to a named repository root directory.
    
    It starts from 'start_path' (usually the script's directory)
    and walks up the filesystem tree until it finds a directory
    matching 'repo_name'.

    Args:
        start_path: The pathlib.Path object to start the search from.
                    Typically the script's own parent directory.
        repo_name: The name of the repository directory to find
                   (e.g., "IncSuperbubble").

    Returns:
        A pathlib.Path object of the found repository root.

    Raises:
        FileNotFoundError: If the repository root directory is not found
                           in any parent path.
    """
    # Start from the resolved, absolute path
    current_dir = start_path.resolve()

    # Loop until we hit the filesystem root (where parent == self)
    while current_dir != current_dir.parent:
        if current_dir.name == repo_name:
            return current_dir  # Found it!
        current_dir = current_dir.parent

    # Final check in case the root directory itself is the repo_name
    if current_dir.name == repo_name:
        return current_dir

    # If we get here, we've hit the filesystem root and found nothing
    raise FileNotFoundError(
        f"Could not find repository root '{repo_name}' in any parent path of {start_path}"
    )

def main():
    """
    Main function to run the test.
    """
    print("Testing repository path discovery...")
    REPO_TO_FIND = "IncSuperbubble"
    
    try:
        # Get the directory containing this test script
        try:
            script_start_path = pathlib.Path(__file__).resolve().parent
        except NameError:
            # Fallback for interactive environments
            script_start_path = pathlib.Path.cwd()
            
        print(f"ℹ️ Starting search from: {script_start_path}")

        # Call the function
        root = find_repo_root(script_start_path, REPO_TO_FIND)
        
        # --- Test Results ---
        print("\n" + "="*30)
        print("✅ SUCCESS: Path discovery test passed.")
        print(f"   Repository Root Found: {root}")
        print(f"   Target Preprocessing Dir: {root / 'Preprocessing'}")
        print("="*30)
        
    except FileNotFoundError as e:
        print("\n" + "!"*30)
        print(f"❌ FAILURE: {e}")
        print("!"*30)
        sys.exit(1)

if __name__ == "__main__":
    main()