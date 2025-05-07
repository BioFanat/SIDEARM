import pandas as pd
import glob
import os
from typing import List, Optional

class CellDataCombiner:
   def __init__(self, directory: str):
       self.directory = directory
       self.combined_df: Optional[pd.DataFrame] = None

   def get_cell_files(self) -> List[str]:
       pattern = os.path.join(self.directory, '*_naked.csv')
       return glob.glob(pattern)

   def process_file(self, filepath: str) -> pd.DataFrame:
       df = pd.read_csv(filepath)
       name = os.path.basename(filepath).replace('_naked.csv', '')
       df['name'] = name
       return df

   def combine_files(self) -> None:
        dfs = [self.process_file(f) for f in self.get_cell_files()]
        self.combined_df = pd.concat(dfs, ignore_index=True)
        self.combined_df.sort_values(['name', 'strand', 'pos'], inplace=True)
        self.combined_df = self.combined_df[['name'] + [col for col in self.combined_df.columns if col != 'name']]
        self.combined_df = self.combined_df[self.combined_df['pos'] != 20.5]

   def save_combined(self, output_path: str) -> None:
       if self.combined_df is None:
           raise ValueError("No data combined yet. Run combine_files() first.")
       self.combined_df.to_csv(output_path, index=False)

def main():
   combiner = CellDataCombiner("/usr/xtmp/bc301/sim_uv_cpd_results")
   combiner.combine_files()
   combiner.save_combined("cpd_naked_collated.csv")

if __name__ == "__main__":
   main()