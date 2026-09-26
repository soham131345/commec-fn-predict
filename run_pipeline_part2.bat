@echo off
REM Copyright (c) 2026 Soham Bhole. All rights reserved.
REM
REM run_pipeline_part2.bat - RESUME script. Use this when Step 3 (data
REM download) already finished but Step 4 got stuck. It shrinks the
REM oversized negatives file, then runs Steps 4-9.
REM
REM Put this .bat in the same folder as run_pipeline.bat (the one with the
REM "scripts" folder), then double-click it.

setlocal
cd /d %~dp0

where conda >nul 2>nul
if not errorlevel 1 goto conda_ok
for %%P in ("%USERPROFILE%\miniconda3" "%USERPROFILE%\Miniconda3" "%USERPROFILE%\anaconda3" "%USERPROFILE%\Anaconda3" "C:\ProgramData\miniconda3" "C:\ProgramData\Miniconda3" "C:\Miniconda3") do (
    if exist "%%~P\Scripts\conda.exe" (
        call "%%~P\Scripts\activate.bat" "%%~P"
        goto conda_ok
    )
)
echo ERROR: conda not found. Run this from Anaconda Prompt instead.
pause
exit /b 1

:conda_ok
call conda activate fnscreen
if errorlevel 1 goto fail

if not exist data\positives.fasta (
    echo ERROR: data\positives.fasta not found - Step 3 never finished.
    echo Run run_pipeline.bat instead.
    pause
    exit /b 1
)

echo ============================================
echo  STEP 3b: shrinking negatives to match positives
echo ============================================
python scripts\subsample_fasta.py data\negatives.fasta data\negatives_balanced.fasta 12439
if errorlevel 1 goto fail

echo ============================================
echo  STEP 4: homology-held-out splits (~15-30 min)
echo ============================================
python scripts\make_homology_splits.py --positives data\positives.fasta --negatives data\negatives_balanced.fasta -o splits --python-fallback
if errorlevel 1 goto fail

echo ============================================
echo  STEP 5: training k-mer baseline
echo ============================================
python scripts\train_kmer_baseline.py --train-fasta splits\train.fasta --train-labels splits\train_labels.csv --val-fasta splits\val.fasta --val-labels splits\val_labels.csv -o models\kmer
if errorlevel 1 goto fail

echo ============================================
echo  STEP 6: training ESM-2 model (slow on CPU, hours)
echo ============================================
python scripts\train_esm_classifier.py --train-fasta splits\train.fasta --train-labels splits\train_labels.csv --val-fasta splits\val.fasta --val-labels splits\val_labels.csv --test-fasta splits\test.fasta --test-labels splits\test_labels.csv -o models\esm
if errorlevel 1 goto fail

echo ============================================
echo  STEP 7: benchmark table
echo ============================================
python scripts\evaluate.py --test-fasta splits\test.fasta --test-labels splits\test_labels.csv --model kmer=models\kmer --model esm=models\esm --esm-test-embeddings models\esm\embeddings_test.npy --esm-test-ids models\esm\ids_test.txt -o benchmark_report.md
if errorlevel 1 goto fail

echo ============================================
echo  STEP 8: drift curves
echo ============================================
python scripts\drift_experiment.py --model-dir models\esm --test-fasta splits\test.fasta --test-labels splits\test_labels.csv -o drift_esm
python scripts\drift_experiment.py --model-dir models\kmer --test-fasta splits\test.fasta --test-labels splits\test_labels.csv -o drift_kmer
if errorlevel 1 goto fail

echo ============================================
echo  STEP 9: demo sequences
echo ============================================
python scripts\make_demo.py --from-test splits\test.fasta splits\test_labels.csv -o demo
if errorlevel 1 goto fail

echo ============================================
echo  ALL DONE. Check benchmark_report.md and drift_esm\drift_curve.png
echo ============================================
pause
exit /b 0

:fail
echo SOMETHING FAILED. Copy the error text above and paste it back for help.
pause
exit /b 1
