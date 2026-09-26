@echo off
REM Copyright (c) 2026 Soham Bhole. All rights reserved.
REM
REM run_pipeline.bat - runs the whole function-prediction pipeline on Windows.
REM
REM BEFORE YOU RUN THIS:
REM   1. Install Miniconda: https://docs.conda.io/en/latest/miniconda.html
REM   2. Put this .bat file in the folder that contains the "scripts" folder.
REM   3. Double-click it (or run it from Anaconda Prompt).
REM      First run takes ~30-60 min (mostly downloads).
REM      ESM training on CPU can take a few hours; on a GPU it's faster.

setlocal
cd /d %~dp0

echo ============================================
echo  STEP 0: checking conda
echo ============================================

where conda >nul 2>nul
if not errorlevel 1 goto conda_ok

for %%P in ("%USERPROFILE%\miniconda3" "%USERPROFILE%\Miniconda3" "%USERPROFILE%\anaconda3" "%USERPROFILE%\Anaconda3" "C:\ProgramData\miniconda3" "C:\ProgramData\Miniconda3" "C:\Miniconda3") do (
    if exist "%%~P\Scripts\conda.exe" (
        call "%%~P\Scripts\activate.bat" "%%~P"
        goto conda_ok
    )
)

echo ERROR: conda not found.
echo Either install Miniconda first, or run this file from "Anaconda Prompt"
echo (find it in the Start Menu).
pause
exit /b 1

:conda_ok
echo conda found.

echo ============================================
echo  STEP 1: creating environment (fnscreen)
echo ============================================
call conda create -y -n fnscreen python=3.11
call conda activate fnscreen
if errorlevel 1 goto fail

echo ============================================
echo  STEP 2: installing packages
echo ============================================
call pip install requests scikit-learn scipy joblib numpy matplotlib torch transformers edlib
if errorlevel 1 goto fail

echo ============================================
echo  STEP 3: downloading training data (UniProt)
echo ============================================
python scripts\fetch_training_data.py -o data --max-positives 10000 --max-negatives 10000
if errorlevel 1 goto fail

echo ============================================
echo  STEP 4: homology-held-out splits
echo ============================================
python scripts\make_homology_splits.py --positives data\positives.fasta --negatives data\negatives.fasta -o splits --python-fallback
if errorlevel 1 goto fail

echo ============================================
echo  STEP 5: training k-mer baseline
echo ============================================
python scripts\train_kmer_baseline.py --train-fasta splits\train.fasta --train-labels splits\train_labels.csv --val-fasta splits\val.fasta --val-labels splits\val_labels.csv -o models\kmer
if errorlevel 1 goto fail

echo ============================================
echo  STEP 6: training ESM-2 model (slow on CPU)
echo ============================================
python scripts\train_esm_classifier.py --train-fasta splits\train.fasta --train-labels splits\train_labels.csv --val-fasta splits\val.fasta --val-labels splits\val_labels.csv --test-fasta splits\test.fasta --test-labels splits\test_labels.csv -o models\esm
if errorlevel 1 goto fail

echo ============================================
echo  STEP 7: benchmark table
echo ============================================
python scripts\evaluate.py --test-fasta splits\test.fasta --test-labels splits\test_labels.csv --model kmer=models\kmer --model esm=models\esm --esm-test-embeddings models\esm\embeddings_test.npy --esm-test-ids models\esm\ids_test.txt -o benchmark_report.md
if errorlevel 1 goto fail

echo ============================================
echo  STEP 8: killer figure (drift curve)
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
echo  ALL DONE.
echo  Results:
echo    benchmark_report.md         <- the model comparison table
echo    drift_esm\drift_curve.png   <- the killer figure
echo    demo\                       <- sequences for the commec demo
echo ============================================
pause
exit /b 0

:fail
echo SOMETHING FAILED. Scroll up, copy the error text, and paste it back for help.
pause
exit /b 1
