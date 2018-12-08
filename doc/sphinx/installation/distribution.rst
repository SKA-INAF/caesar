============
Distribution
============

To support users having problems to install CAESAR from source in their system we provide Docker and Singularity container recipe files and pre-built images with all software installed.
 

---------------------
Docker base container
---------------------

We assume here that Docker is available in your system. If not, install it following the instruction at https://docs.docker.com.

Recipe files to build a base Docker image with all CAESAR dependencies installed can be downloaded at: https://github.com/SKA-INAF/caesar-base-docker.
For example the following command builds a Docker image on your system with name `caesar/base` (tag=latest) from `Dockerfile.xenial`: 


``docker build -t 'caesar/base:latest' -f 'Dockerfile.xenial' .``


A pre-built Docker image (size approximately ~16 GB), ready to be used, can be download from `Docker Hub` using the following command:


``docker pull sriggi/caesar-base``


The image is built over a Ubuntu 16:04 (xenial) base OS image. Software dependencies are installed under the ``/opt/Software`` directory.

A Docker container with CAESAR installed will be provided in the future.


---------------------
Singularity container
---------------------

We assume here that Singularity is available in your system. If not, install it from https://singularity.lbl.gov.

Recipe files to build the CAESAR singularity image can be downloaded at: https://github.com/SKA-INAF/caesar-singularity. Singularity image is built using the base Docker image
described above.
For example the following command builds a production container (e.g. non writable) on your system with name `caesar.simg` from a recipe file `Singularity.xenial`:


``sudo singularity build caesar.simg Singularity.xenial``


A pre-built container image (size approximately ~2 GB), ready to be used, can be download from `Singularity Hub` using the following command:


``singularity pull --name caesar.simg shub://SKA-INAF/caesar-singularity:xenial``


To enter the container type:


``singularity shell caesar.simg``


CAESAR is installed in the ``/opt/Software/CAESAR/trunk/`` container directory. Software dependencies are installed under the ``/opt/Software`` directory. 
CAESAR applications are available for batch processing as Singularity apps. To list available apps type:


``singularity apps caesar.simg``


To run one of these applications (say with name `APP_NAME` and with input arguments `APP_ARGS`) type:


``singularity run --app APP_NAME caesar.simg APP_ARGS``


For example source finding (see Tutorial section) can be run inside container with:


``singularity run --app sfinder caesar.simg --config=[FILE]``


