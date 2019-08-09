"""Controls the Monte Carlo runs."""
import os
import sys
import numpy as np
from astropy.time import Time
from .Output import Output


def delete_files(filelist):
    """Delete output files and remove them from the database.
    
    **Parameters**
    
    filelist
        List of files to remove. This can be found with Inputs.findpackets()
        
    **Returns**
    
    No outputs.
    
    """
    from .database_connect import database_connect

    with database_connect() as con:
        cur = con.cursor()

        for f in filelist:
            # Delete the file
            print(f)
            if os.path.exists(f):
                os.remove(f)

            # Remove from database
            cur.execute('''SELECT idnum FROM outputfile
                           WHERE filename = %s''', (f, ))
            idnum = cur.fetchone()[0]

            cur.execute('''DELETE FROM outputfile
                           WHERE idnum = %s''', (idnum, ))
            cur.execute('''DELETE FROM geometry
                           WHERE geo_idnum = %s''', (idnum, ))
            cur.execute('''DELETE FROM sticking_info
                           WHERE st_idnum = %s''', (idnum, ))
            cur.execute('''DELETE FROM forces
                           WHERE f_idnum = %s''', (idnum, ))
            cur.execute('''DELETE FROM spatialdist
                           WHERE spat_idnum = %s''', (idnum, ))
            cur.execute('''DELETE FROM speeddist
                           WHERE spd_idnum = %s''', (idnum, ))
            cur.execute('''DELETE FROM angulardist
                           WHERE ang_idnum = %s''', (idnum, ))
            cur.execute('''DELETE FROM options
                       WHERE opt_idnum = %s''', (idnum, ))
            print(f'Removed {idnum}: {os.path.basename(f)} from database')

            cur.execute('''SELECT idnum, filename FROM modelimages
                           WHERE out_idnum = %s''', (idnum, ))
            for mid, mfile in cur.fetchall():
                cur.execute('''DELETE from modelimages
                               WHERE idnum = %s''', (mid, ))
                if os.path.exists(mfile):
                    os.remove(mfile)

            cur.execute('''SELECT idnum, filename FROM uvvsmodels
                           WHERE out_idnum = %s''', (idnum, ))
            for mid, mfile in cur.fetchall():
                cur.execute('''DELETE from uvvsmodels
                               WHERE idnum = %s''', (mid, ))
                if os.path.exists(mfile):
                    os.remove(mfile)
