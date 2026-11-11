# Changes.docx implementation

1. Retained `equal_total_250`, including original source CSVs and canonical Figure 1/2 filenames.
2. Added `equal_exposure_c5`: endpoints 200, 300, 400, 500, 600 for N=40,60,80,100,120.
3. Added absolute and ancestor-subtracted breadth-versus-exposure figures, faceted by mutation kernel and chi, comparing all five redundancies through x=5.
4. Both endpoint tables use identical metric functions, and both endpoints come from the same seeded history. No endpoint is extrapolated.
5. The primary execution targets five seeds per cell. A 20-seed configuration and deterministic extension registry are implemented; expanded execution status is recorded in the execution summary. Original five seeds are retained, and validated completed runs are skipped.

See `Changes.transcribed.txt` for the full co-author wording. All scientific differences are stopping/measurement analyses; no new biological model condition was introduced.
