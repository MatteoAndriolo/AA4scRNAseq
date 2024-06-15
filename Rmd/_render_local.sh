#!/bin/bash

# docker exec agitated_nobel /app/Rmd/_main.sh Melanoma
docker run --rm -v $(pwd):/app -w /app myrocker:latest /app/Rmd/_main.sh Melanoma