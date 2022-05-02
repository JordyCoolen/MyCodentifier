# ${1} <filename> example: mycoprofiler
# ${2} <imagename> example: mycoprofiler:version

docker run -v /var/run/docker.sock:/var/run/docker.sock \
  --privileged \
  -it \
  --rm \
  -v /Users/jordycoolen/PycharmProjects:/output \
  quay.io/singularity/docker2singularity \
  --name ${1} \
  ${2} /bin/bash