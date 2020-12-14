docker build -t jmonlong-svnicu .

docker tag jmonlong-svnicu quay.io/jmonlong/svnicu:0.3
docker push quay.io/jmonlong/svnicu:0.3

docker run -it -v `pwd`:/app -w /app -u `id -u $USER` jmonlong-svnicu
