#!/bin/bash
squeue -u andrioloma -h -o %i | xargs scancel