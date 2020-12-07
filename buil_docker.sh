docker build -t jmonlong-svnicu .

docker tag jmonlong-svnicu quay.io/jmonlong/svnicu
docker push quay.io/jmonlong/svnicu

docker run -it -v `pwd`:/app -w /app -u `id -u $USER` jmonlong-svnicu
