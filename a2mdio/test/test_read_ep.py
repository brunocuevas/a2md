from a2mdio.qm import GaussianLog

gl = GaussianLog(
    file='test_08.g09', verbose=True
)
gl.read()
ep = gl.seek_ep()

# for i, v in enumerate(ep):
#     print("{:06d} {:8.4e}".format(i, v))
