#!/bin/bash
src_dir=$(dirname "$0")

#download Gene locations, build 37
wget https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq/download -O ${src_dir}/NCBI37.3.zip
gunzip ${src_dir}/NCBI37.3.zip

#donwload B files
wget https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download -O ${src_dir}/g1000_eur.zip 
gunzip ${src_dir}/g1000_eur.zip 
