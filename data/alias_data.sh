# decompress all .gz
alias dcgzip='find . -type f -name "*.gp" -exec gzip -dk {}'
alias untar='tar -czvf' # insert later output file and files to compress
alias maketar='tar -xzvf' #output .tar.gz and files to compress