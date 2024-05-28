#!/bin/bash

# sync local out/ folders with folders from server

rsync -avz --exclude '.*' --exclude '*/.*' --exclude '*.rds' andrioloma@login.dei.unipd.it:/home/andrioloma/MasterThesis/out .

rsync -avz --include="*/" --include='**/*.html' --include='*_files/**/*' --exclude='*' andrioloma@login.dei.unipd.it:/home/andrioloma/MasterThesis/Rmd .