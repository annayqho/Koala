""" Log onto Kowalski """

from penquins import Kowalski

def logon():
    """ Log onto the database """
    username = 'ah'
    password = 'TetraodonInsists'
    s = Kowalski(
            protocol='https', host='kowalski.caltech.edu', port=443,
            verbose=False, username=username, password=password)
    print("logged on")
    return s
