## Instructions to run the test case

1. In `preprocessing` directory, run `process_tags.m`. It reads raw tag data file downloaded from the tag, located in `raw_tag_data` directory, and outputs a .mat file (e.g., `7_raw.mat`).
  * Specify in `ptags` array the tag numbers you want to preprocess.
2. Run `run_tag.m` in this directory.
  * Specify in `ptags` array the tag numbers you want to run.
  * `do_parts(1)` generates the 2-D grid required for HMM computation and is only needed when running the first tag from the batch.
