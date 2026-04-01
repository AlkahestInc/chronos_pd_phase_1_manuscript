import random
import math

def get_assay_ids(assay_ids_file):
  """
  Returns assay_ids from file.
  
  Args:
    assay_ids_file (str): Path to file containing Assay IDs (one per line)
   
  Returns:
    list[str]: List of assay_ids.
   
  Raises:
    FileNotFoundError: If file doesn't exist (with helpful message)
  """
  try:
    with open(assay_ids_file) as f:
        assay_ids = [line.strip() for line in f if line.strip()]
        if not assay_ids:
          raise ValueError(f"File {assay_ids_file} is empty")
        return assay_ids
  except FileNotFoundError:
    # More explicit error for debugging
    raise FileNotFoundError(
      f"Assay IDs file not found: {assay_ids_file}\n"
      f"Run 'snakemake extract_assay_ids' first, or use a checkpoint."
    )

def select_equally_spaced(items, n):
  """
  Select n equally spaced items from the list.
  """
  if n is None or n >= len(items):
    return items
  step = len(items) / n
  return [items[math.floor(i * step)] for i in range(n)]
