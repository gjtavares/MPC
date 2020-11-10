function [u_fut]=alg_controle(Ka,w_ref,fu,D,p,fd)

u_fut=Ka*(w_ref-fu-D*p-fd);