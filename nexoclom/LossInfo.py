import numpy as np
import psycopg2
import astropy.units as u
import astropy.units.astrophys as ua

class LossInfo():
    def __init__(self, atom, lifetime, aplanet, database):
        # Initialization
        self.photo = 0.
        self.eimp = 0.
        self.chX = 0.
        self.reactions = []

        if lifetime.value < 0:
            self.photo = np.abs(1./lifetime.value)
            self.reactions = 'Generic photo reaction.'
        else:
            con = psycopg2.connect(host='localhost', database=database)
            cur = con.cursor()

            # Photo rate adjusted to proper heliocentric distance
            cur.execute('''SELECT reaction, kappa
                           FROM photorates
                           WHERE species=%s and bestvalue=True''',
                        (atom, ))
            if cur.rowcount == 0:
                print('No photoreactions found')
            else:
                rows = cur.fetchall()
                for r in rows:
                    self.reactions.append(r[0])
                    self.photo += r[1]/aplanet**2

            # Electron impact

            # Charge exchange

        if len(self.reactions) == 0:
            self.reactions = None

    def __len__(self):
        return len(self.reactions) if self.reactions is not None else 0

    def __str__(self):
        if len(self) == 0:
            print('No reactions included')
        elif len(self) == 1:
            print('Included Reaction: {}'.format(self.reactions[0]))
        else:
            print('\tIncluded Reactions: {}'.format(tuple(self.reactions)))
        if self.photo != 0:
            print('Photo Rate = {:0.2e} s'.
                  format(self.photo))
        if self.eimp != 0:
            print('Electron Impact Rate = {:0.2e} UNIT'.
                  format(self.eimp))
        if self.chX!= 0:
            print('Charge Exchange Rate = {:0.2e} UNIT'.
                  format(self.chX.value))

        return ''
