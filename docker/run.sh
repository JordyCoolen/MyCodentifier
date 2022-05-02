# ${1} <containername> example: mycodentifier
# ${2} <imagename> example: mycodentifier:1.0

docker run \
  -it \
  --rm \
  --name ${1} \
  --mount type=bind,source=$PWD,target=/workflow \
  ${2} /bin/bash