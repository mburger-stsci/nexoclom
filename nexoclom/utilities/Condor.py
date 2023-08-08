import os
import sys
import time
import random
try:
    import htcondor
except:
    pass


class Condor:
    def __init__(self, njobs=1000, delay=5, nice=None):
        self.njobs = njobs
        self.delay = delay
        self.nice = nice
        self.schedd = htcondor.Schedd()

    def nrunning(self):
        return len(list(self.schedd.xquery()))
    
    def submit(self, command, **kwargs):
        ct = 0
        # If more than njobs are running, wait until some complete
    
        while self.nrunning() >= self.njobs:
            print('Waiting 30 s')
            time.sleep(30)
        
        time.sleep(random.random()*self.delay)
        
        logfile = kwargs.get('logfile', 'logs/log.log')
        outlogfile = kwargs.get('outlogfile', 'logs/out.out')
        errlogfile = kwargs.get('errlogfile', 'logs/err.err')
        if not os.path.exists(os.path.dirname(logfile)):
            os.makedirs(os.path.dirname(logfile))
    
        submit = {'universe': 'vanilla',
                  'executable': sys.executable,
                  'getenv': '*',
                  'log': logfile,
                  'output': outlogfile,
                  'error': errlogfile}

        if 'arguments' in kwargs:
            submit['arguments'] = kwargs['arguments']
        else:
            pass

        if 'request_memory' in kwargs:
            submit['request_memory'] = kwargs['request_memory']
        else:
            pass

        requirements = []
        if 'machine' in kwargs:
            requirements.append(
                f'''TARGET.Machine == "{kwargs['machine']}.stsci.edu" ''')
        else:
            pass

        if len(requirements) > 0:
            submit['requirements'] = '&&'.join(requirements)
        else:
            pass
            
        job = htcondor.Submit(submit)
        cluster = self.schedd.submit(job).cluster()
    
        # with schedd.transaction() as txn:
        #     print(f'Submitting job {sub.queue(txn)}')
        #     jobnumber = sub.queue(txn)
      
        if self.nice is not None:
            os.system(f'renice {self.nice} -u mburger')
            
        return cluster

    def nCPUs(self):
        collector = htcondor.Collector()
        query = collector.query()
        num, names, names_str = 0, [], []
        for ad in query:
            if (('DetectedCpus' in ad.keys()) and (ad.lookup('Name') not in names) and
                ('@' not in ad.lookup('Name'))):
                names.append(ad.lookup('Name'))
                names_str.append(ad.lookup('Name').eval().split('.')[0])
                num += int(ad.lookup('DetectedCpus'))
            else:
                pass
        
        return num, names_str

    def activeJobs(self):
        return set(int(ad.lookup('ClusterId')) for ad in htcondor.Schedd().xquery())

    def n_to_go(self, cluster_list):
        active = self.activeJobs()
        return len(cluster_list) - len(set(cluster_list) - active)

    def machines_in_use(self):
        hosts = set(str(ad.lookup('RemoteHost'))
                    for ad in htcondor.Schedd().xquery()
                    if 'RemoteHost' in ad.keys())
        machines = set(x.split('@')[1].split('.')[0] for x in hosts)
        return machines
