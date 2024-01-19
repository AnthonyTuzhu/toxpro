# Welcome to ToxiVerse Project

![repo version](https://img.shields.io/badge/Version-v.%200.2.1-green)

*Added regression QSAR modeling module (experimental)
Fixed some bugs
### Running in PyCharm (IDE)


To open this project on your computer, you need to set up the conda env by using the `requirements.txt` file to install
the necessary python dependencies.

First, create a new conda env

```
conda env create -n toxpro
```

then you need to activate

```
conda activate toxpro
```

close the terminal in IDE (e.g., Pycharm), open any .py file under toxpro/app and re-open the terminal in IDE

then use pip to install dependencies, for example:

```
...\anaconda3\envs\toxpro\python.exe -m pip install -r requirements.txt
```

then download the database folders on Google Drive because Github lacks support to sharing large files:
https://drive.google.com/drive/folders/1gBxySt2uk6XNflRyA7mf5y-zqY0p2ErS?usp=drive_link (data) and
https://drive.google.com/drive/folders/1a-lKJ-JwjyY-v3J4DpRgCP7prCEf2WQP?usp=drive_link (instance).

Since google drive download folders as zipped file, drag the folders to the work respiratory location (e.g., toxpro, replace existing ones if there are).

Remember to set up class Config in the config.py:
```
MASTER_DB_FILE = os.getenv('MASTER_DB_FILE')
BIOASSAYS = os.getenv('BIOASSAYS')
BIOPROFILE_DIR = os.getenv('BIOPROFILE')
```

Intiate the Flask app by __init__.py


Getting started using the Flask tutorial: https://flask.palletsprojects.com/en/2.0.x/tutorial

https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xix-deployment-on-docker-containers

### CAUTION : Docker is ESSENTIAL for QSAR module running: https://www.docker.com/blog/how-to-use-the-redis-docker-official-image/

After installed the Docker app, you will need to run toxpro/docker_compose.yml to compose the environment and then click on "5000:5000" port in the "container" session of Docker

Here is a tutorial for using it with flask: https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xxii-background-jobs

For master-database, since bioprofile is not completed, some table may raise ajax error. But it won't affect other information.

## Run app in command lines

If you want to run the app in command line, first cd the location of this app, for example: 
```
C:\Users\Anthony>cd toxpro
```

Similarily, install requirements.txt
```
python -m pip install -r requirements.txt
```

Download "data" and "instance" folder and put them under the app folder, set up class Config in the config.py, and run the flask app
```
flask run
```

To run Docker, use the docker-compose.yml:

```
docker-compose.exe -f docker-compose.yml -p toxpro up -d
```

## Redis Queue, Docker and GUNICORN


Redis Queue is a way to handle background tasks/jobs.  It requires a seperate server and installation outside
of python.



There are two commands that need to be run to get the redis server runnign using docker. The first is to start the redis
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

push all the code with somethign like this:

```commandline
scp -r ~/projects/toxpro/* toxpro@192.241.131.84:/home/toxpro
```

Log in and using docker compose to build containers

```commandline
docker compose up -d
```


## Docker commands

remove all containers
```dockerfile
sudo docker system prune -a
```

force recreate of building image and containers
```dockerfile
docker compose -f docker-compose-do.yml up
```
or this
```dockerfile
docker compose -f docker-compose-do.yml up -d --build