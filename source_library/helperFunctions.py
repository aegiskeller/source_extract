import astropy.io.fits as fits
import sqlite3


# create a function that prints the contents of a fits table
def printTable(filename):
    # open the table file
    hdul = fits.open(filename)
    # print the header
    print(hdul[0].header)
    # print the data
    print(hdul[1].data)
    # close the table file
    hdul.close()


# create or connect to a SQLite database
def connect_to_db(db_name):
    """Connect to the SQLite database and return the connection object."""
    conn = sqlite3.connect(db_name)
    return conn


def create_table(conn, create_table_sql):
    """Create a table from the create_table_sql statement."""
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Exception as e:
        print(e)


# a function to create a table 'wcs' in the database
def create_wcs_table(conn):
    """Create a table 'wcs' in the database.
    ra_center, dec_center, pixscale, orientation, fieldw, fieldh, fieldunits
    """
    sql = """CREATE TABLE IF NOT EXISTS wcs (
          id integer PRIMARY KEY,
          ra_center real NOT NULL,
          dec_center real NOT NULL,
          pixscale real NOT NULL,
          orientation real NOT NULL,
          fieldw real NOT NULL,
          fieldh real NOT NULL,
          fieldunits text NOT NULL
    );"""
    create_table(conn, sql)


def resolve_coordinates(target):
    """Resolve the coordinates of a target using the SIMBAD service."""
    from astroquery.simbad import Simbad

    result_table = Simbad.query_object(target)
    # sometimes the result_table has 'RA' sometimes 'ra'
    if "RA" not in result_table.colnames:
        ra = result_table["ra"][0]
        dec = result_table["dec"][0]
    else:
        ra = result_table["RA"][0]
        dec = result_table["DEC"][0]
    return ra, dec
