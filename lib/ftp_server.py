import subprocess
import argparse
from pyftpdlib.authorizers import DummyAuthorizer
from pyftpdlib.handlers import FTPHandler
from pyftpdlib.servers import FTPServer




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='show mol file')
    parser.add_argument('--username', "-U", type=str, default="root", help='username default: root')
    parser.add_argument('--password', "-P", type=str, default="abc123", help='password')
    parser.add_argument('--path', type=str, default=r'C:\Users\admin', help='default path')
    parser.add_argument('--port', "-p", type=int, default=2122, help='port')
    args = parser.parse_args()
    std = subprocess.Popen("ipconfig", shell=True)
    std.wait(timeout=1)
    print(f"uname: {args.username}; pwd: {args.password}; port: {args.port}; path: {args.path}")

    #实例化虚拟用户，这是FTP验证首要条件
    authorizer = DummyAuthorizer()
    # python -m pyftpdlib -p 2122 -u root -P "********"
    #添加用户权限和路径，括号内的参数是(用户名， 密码， 用户目录， 权限)
    authorizer.add_user(args.username, args.password, args.path, perm='elradfmw')

    #添加匿名用户 只需要路径
    #初始化ftp句柄
    handler = FTPHandler
    handler.authorizer = authorizer

    #监听ip 和 端口
    server = FTPServer(('0.0.0.0', args.port), handler)

    #开始服务
    server.serve_forever()