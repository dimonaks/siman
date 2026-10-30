# SSHTools module
#
# Andriy Zhugayevych (azh@ukr.net), Sergei Matveev(matseralex@yandex.ru)
# www.zhugayevych.me/python/SSHTools/index.htm
# created 21.08.2014, modified 12.09.2016

import paramiko
import os
import sys
import hashlib
import shutil

class SSHTools:
	user=""; host=""; pkeypath=""; pkey=""; port=22
###############################################################################################
	def setup( self, user="", host="", pkey="", port=22 ):
		if user:
			self.user = user
		if host:
			self.host = host
		if pkey:
			self.pkeypath = pkey
		self.port = port
		self.ssh = paramiko.SSHClient()
		self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		self.pkey = paramiko.RSAKey.from_private_key_file(pkey)
###############################################################################################
	def run( self, command, noerror=False, printout=False ):
		if printout:
			print("-- SSHTools(): command:",command)
		# print("command:",command)
		self.ssh.connect(self.host,username=self.user,pkey=self.pkey,port=self.port)
		stdin,stdout,stderr = self.ssh.exec_command(command)
		out=stdout.readlines()
		err=stderr.readlines()
		self.ssh.close()
		if err:
			if noerror:
				if printout:
					print("-- SSHTools(): stderr:","".join(err))
				return ("".join(err)).strip()
			else:
				raise Warning("-- SSHTools(): stderr: "+"".join(err))
		
		if printout:
			if len(out)==1:
				print("-- SSHTools(): output:",out[0].rstrip())
			else:
				print("-- SSHTools(): output:")
				for s in out:
					print(s.rstrip())
		return ("".join(out)).strip()
###############################################################################################
	def put( self, source, dest ):
		self.ssh.connect(self.host,username=self.user,pkey=self.pkey,port=self.port)
		sftp = self.ssh.open_sftp()
		
		# print(source, dest)
		# print(os.getcwd())
		# print(os.path.exists(source) )
		sftp.put(source,dest)

		sftp.close()
		self.ssh.close()
		return
###############################################################################################
	def _get_file(self, sftp, source, dest, progress=False):
		# Sequential copy. Paramiko's pipelined sftp.get stalls or stops
		# on files larger than about 200 MB when the SSH window is small.
		try:
			channel = sftp.get_channel()
			if channel is not None and getattr(channel, 'lock', None) is not None:
				channel.lock.acquire()
				try:
					channel.in_window_size = max(getattr(channel, 'in_window_size', 0), 2147483647)
					channel.out_window_size = max(getattr(channel, 'out_window_size', 0), 2147483647)
				finally:
					channel.lock.release()
		except Exception:
			pass
		try:
			size = sftp.stat(source).st_size
		except Exception:
			size = 0
		name = os.path.basename(source)
		copied = 0
		shown = -1
		with sftp.open(source, 'rb') as remote_file, open(dest, 'wb') as local_file:
			remote_file.MAX_REQUEST_SIZE = 1024 * 1024
			while True:
				block = remote_file.read(8 * 1024 * 1024)
				if not block:
					break
				local_file.write(block)
				copied += len(block)
				if not progress:
					continue
				percent = int(100 * copied / size) if size else 100
				if percent != shown:
					shown = percent
					filled = 30 * percent // 100
					bar = '#' * filled + '-' * (30 - filled)
					sys.stdout.write('\r%s  [%s] %3d%%' % (name, bar, percent))
					sys.stdout.flush()
		if progress and shown >= 0:
			sys.stdout.write('\n')
			sys.stdout.flush()

	def get( self, source, dest, progress=False ):
		# print(source, dest )
		self.ssh.connect(self.host,username=self.user,pkey=self.pkey,port=self.port)
		sftp = self.ssh.open_sftp()
		self._get_file(sftp, source, dest, progress=progress)
		sftp.close()
		self.ssh.close()
		return
###############################################################################################
	def sget( self, source, dest, keeplocal=False ):
		self.ssh.connect(self.host,username=self.user,pkey=self.pkey,port=self.port)
		sftp = self.ssh.open_sftp()
		self._get_file(sftp, source, dest)
		sftp.close()
		self.ssh.close()
		s=self.run("md5sum "+source)
		md5_remote=str(s).split()[0]
		md5_local=hashlib.md5(open(dest,'rb').read()).hexdigest()
		if md5_local!=md5_remote:
			print("md5 sums are different for",source,"and",dest,": ",md5_remote,"!=",md5_local)
			if not(keeplocal):
				os.remove(dest)
				print("Local copy is removed")
			raise Warning("Checksum error")
		return 0

###############################################################################################
	def fexists(self, filename):
#		return isinstance ( self.run ("ls "+ filename, noerror=True), str)
		return filename == self.run("ls " + filename, noerror = True)
