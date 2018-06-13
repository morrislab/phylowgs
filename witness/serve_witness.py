from SimpleHTTPServer import SimpleHTTPRequestHandler
import SocketServer
import argparse
import StringIO
import gzip
import zipfile

class WitnessServerHandler(SimpleHTTPRequestHandler):
  def do_GET(self):
    #Extend the do_GET function so that if the client requests a summ or muts json file,
    #it loads it, gzip compresses it, and then sends it. In addition, if the client 
    #requests a specific file out of mutass.zip, then extract that file, gzip it, and
    #send it. This will save us a lot of time by not requiring many large files be sent
    #over the network.
    pathToFile = self.path
    if('summ.json' in pathToFile) | ('muts.json' in  pathToFile):
      f = file(pathToFile[1:], 'r')
      self._SendFile(f)
    elif 'mutass.zip' in pathToFile:
      #Extract the file of interest, gzip encode it, and send it.
      fileToExtract = pathToFile.split('/')[-1]
      zipPath = pathToFile[0:len(pathToFile)+len(fileToExtract)]
      mutassZip = zipfile.ZipFile(zipPath, mode='r')
      f = mutassZip.read(fileToExtract)
      self._SendFile(f)
    else:
      SimpleHTTPRequestHandler.do_GET(self)
  
  def _SendFile(self,f):
    fdat = f.read();
    tosend = self._GzipEncode(fdat)
    self.send_response(200)
    self.send_header("Content-type", "text/html")
    self.send_header("Content-length", len(str(tosend)))
    self.send_header("Content-Encoding", "gzip")
    self.end_headers()
    self.wfile.write(tosend)
    self.wfile.flush()

  def _GzipEncode(self, content):
    out = StringIO.StringIO()
    f = gzip.GzipFile(fileobj=out, mode='w', compresslevel=5)
    f.write(content)
    f.close()
    return out.getvalue()

def parse_args():
  parser = argparse.ArgumentParser(
    description='Run a server to run witness on.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument(dest='port', default=8000, type=int,
    help='Which port the server handler should listen to. Default = 8000')
  args = parser.parse_args()
  return args

def main():
  args = parse_args()
  httpd = SocketServer.TCPServer(("", args.port), WitnessServerHandler)
  print("Serving at port " + str(args.port))
  httpd.serve_forever()

if __name__ == "__main__":
    main()