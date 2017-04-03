import sys
import os

#exts = ['.cs','.aspx','.cs','.js','.css','.ashx','.master','.ascx','.svc','.config']
exts = ['.py', '.f']
fortranPath = "C:\\Users\\12442063\\Dropbox\\Share with MAC\\FLiES\\FLiES_v2.43\\src"

def read_line_count(fname):
    count = 0
    for file_line in open(fname).readlines():
        count += 1
    return count

if __name__ == '__main__':

    count = 0
    fcount = 0
    for root,dirs,files in os.walk(fortranPath):# os.getcwd()
        for f in files:
            # Check the sub directorys
            fname = (root + '/'+ f).lower()
            print(fname)
            ext = f[f.rindex('.'):]
            try:
                if(exts.index(ext) >= 0):
                    fcount += 1
                    c = read_line_count(fname)
                    count += c
            except:
                pass

    print('file count:' + str(fcount))
    print('count:' + str(count))
