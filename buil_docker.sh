docker build -t jmonlong-svnicu .

docker tag jmonlong-svnicu quay.io/jmonlong/svnicu:0.5
docker push quay.io/jmonlong/svnicu:0.5

docker run -it -v `pwd`:/app -w /app -u `id -u $USER` jmonlong-svnicu
