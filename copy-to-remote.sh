#!/bin/bash
copy_command="rsync -avz"
copy_args="--exclude=.* --exclude=build --exclude=cmake*  --delete --delete-excluded"

# set default path and server
local_code_path="$(dirname "$(readlink -f "$0")")"
remote_proj_name=VtsReader
remote_code_path="shannon:/lustre/home/huachengli/${remote_proj_name}/"

if [ $1 = "s13" ]; then
	remote_code_path="s13:/home2/huachengli/${remote_proj_name}/"
elif [ $1 = "shannon" ]; then
	remote_code_path="shannon:/lustre/home/huachengli/${remote_proj_name}/"
elif [ $1 = "mars" ]; then
	remote_code_path="mars:/public3/home/spa/salecvtsreader_dev/${remote_proj_name}/"
elif [ $1 = "crater" ]; then
	remote_code_path="crater:/public/home/huachengli/${remote_proj_name}/"
else
	echo "USE DEFAULT SERVER"
fi
echo "COPY ${local_code_path} TO ${remote_code_path}"
$copy_command $copy_args ${local_code_path}/ $remote_code_path
