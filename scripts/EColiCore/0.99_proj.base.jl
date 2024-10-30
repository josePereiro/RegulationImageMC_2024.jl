# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
# proc
_procdir(args...) = projpath(PROJ, ["data", string(PROJVER)], args...)

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
# db
_blobsdir(args...) = _procdir(["blob"], args...)

nothing

