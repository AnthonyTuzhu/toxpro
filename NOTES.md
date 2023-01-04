### BASIC SETUP

To set up the conda env, you need to use the `requirements.txt` file to install
the necessary python dependencies.

First, create a new conda env

```
conda env create -n toxpro
```

then you need to activate

```
conda activate toxpro
```
then use pip to install dependencies

```
pip install -f requirements.txt
```



Getting started using the Flask tutorial: https://flask.palletsprojects.com/en/2.0.x/tutorial

https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xix-deployment-on-docker-containers


command for pushing image to docker hub

docker tag nanodb:latest danrusso/nanodb:latest
docker push danrusso/nanodb


command for starting docker on digital ocean

need to remove image first to make sure it pulls from repo after update

sudo docker run -p 80:5000 -v /c/Users/russod/Projects/toxpro/logs:/home/toxpro/logs -v /c/Users/russod/Projects/toxpro/data:/home/tox/data/ --env NANO_FILE_DIR=/home/tox/data --name toxpro -d -t toxpro

I was getting a exec ./boot.sh: no such file or directory


turns out i was getting removing windows style line endings this solved the problem:
https://unix.stackexchange.com/questions/79702/how-to-test-whether-a-file-uses-crlf-or-lf-without-modifying-it
https://stackoverflow.com/questions/11680815/removing-windows-newlines-on-linux-sed-vs-awk

helpful notes:

SQLALCHEMY and app factories:https://flask.palletsprojects.com/en/2.2.x/patterns/appfactories/

started using flask migrations: https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-iv-database

flask db init to intializaize db

Updates to the `db_models.py` classes need to be reflected in the database.  this requires two commands to be run.


## Redis Queue

Redis Queue is a way to handle background tasks/jobs.  It requires a seperate server and installation outside
of python.

Docker should be good for this: https://www.docker.com/blog/how-to-use-the-redis-docker-official-image/

Here is a tutorial for using it with flask: https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xxii-background-jobs


There are two commands that need to be run to get the redis server runnign using docker.  The first is to start the redis
container using:

```
docker run --name my-redis -p 6379:6379 redis
```

then in a separate terminal need to start the RQ worker.  This needs to be done with the ToxPro env activated and in the
directory root of the project (i.e., same level as `app/`)

```
rq worker toxpro-tasks
```
### Docker, Redis and rq

Docker, redis and rq were not trivial to get working together, but [this](https://github.com/fcakyon/flask-redis-docker) template
structure and [this](https://blog.abbasmj.com/implementing-redis-task-queues-and-deploying-on-docker-compose) blog post were 
pretty helpful.  This biggest challenge was getting recognizing that the container names (e.g., `toxpro-redis-1`) was required
for connecting from the app.  There was also an issue with naming the workers and the queue (apparently these are different).


## Email

was able to set up the email using google and this: https://stackoverflow.com/questions/72478573/how-to-send-an-email-using-python-after-googles-policy-update-on-not-allowing-j


## Docker and GUNICORN

Runing the app in the docker container was very slow unitl I changed the setings in `.\boot.sh` as per this blog:
https://pythonspeed.com/articles/gunicorn-in-docker/.

## Deployment on Digital Ocean

-Created a project with a droplet at the cheapest price.  I used the 1-click docker
installed droplet as outlined [here](https://www.digitalocean.com/community/tutorials/how-to-use-the-docker-1-click-install-on-digitalocean)

I then created images of all the containers and used `docker compose push` to push them
to my repo [here are the docs](https://docs.docker.com/engine/reference/commandline/compose_push/) for that command.

I then copied the `.env` and `docker-compose.yml` files over to the droplet and did the `docker compose pull`. 

I did have to create an instance of the of the `instance` folder to map the sqlite and copy that 
over as well.  

There is probably a better way to do this. 

Well the above did not work, due to the M1 chips and building on a mac vs linux.  The temp solution
now is just to copy all the code to and build on the droplet. 

Need to create a new user toxpro because I was getting errors writing to the SQLite database
due to permissions.  Outlined [here](https://www.digitalocean.com/community/tutorials/how-to-create-a-new-sudo-enabled-user-on-ubuntu-18-04-quickstart)

Then, I needed to intall docker compose as outline [here](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-compose-on-ubuntu-22-04).