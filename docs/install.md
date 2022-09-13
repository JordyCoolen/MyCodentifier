## installation notes:
### OPTION1 (build docker image from scratch)
```bash
install docker https://www.docker.com/products/docker-desktop
git clone https://github.com/JordyCoolen/MyCodentifier.git
cd MyCodentifier
docker build --rm -t jonovox/mycodentifier:1.0 ./
```

### or

### OPTION2 (pull docker image from docker hub)
```bash
install docker https://www.docker.com/products/docker-desktop
git clone https://github.com/JordyCoolen/MyCodentifier.git
cd MyCodentifier
docker pull jonovox/mycodentifier:1.0
```

### CONDA environments (prebuild)
download the newest release of the pipeline via:
https://surfdrive.surf.nl/files/index.php/s/azkU2t09zY1r9qB
and extract into the Mycodentifier folder

## Databases
- currently the databases are not provided
- we will provide in future versions the databases
- please contact us for more information

## create singularity .sif out of docker image
```bash
bash docker/create_singularity_image.sh <filename> <imagename>
```
example:
```bash
bash docker/create_singularity_image.sh mycodentifier_1_0 mycodentifier:1.0
```