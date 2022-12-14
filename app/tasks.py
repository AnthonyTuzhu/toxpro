# this module is outlined
# here: https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xxii-background-jobs
from app import create_app
import sys, time

app = create_app()
app.app_context().push()

def build_qsar_models():
    try:
        pass
    except:
        app.logger.error('Unhandled exception', exc_info=sys.exc_info())
    #finally:

def example(seconds):
    print('Starting task')
    for i in range(seconds):
        print(i)
        time.sleep(1)
    print('Task completed')