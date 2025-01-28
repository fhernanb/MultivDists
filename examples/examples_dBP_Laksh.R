# Pruebas para chequear vectorizacion de la funcion de densidad
dBP_Laksh(c(0, 0), l1=3, l2=4, alpha=1)
dBP_Laksh(c(1, 0), l1=3, l2=4, alpha=1)
dBP_Laksh(c(0, 1), l1=3, l2=4, alpha=1)

muestra <- matrix(c(0, 0,
                    1, 0,
                    0, 1), ncol=2, byrow=TRUE)
muestra
dBP_Laksh(x=muestra, l1=3, l2=4, alpha=1)

