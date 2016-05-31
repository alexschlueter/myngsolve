objects = myHOElement.o myHOFESpace.o

%.o : %.cpp
	${NETGENDIR}/ngscxx -I. -c $? -o $@


libmyngsolve.so : $(objects)
	${NETGENDIR}/ngsld -shared $(objects) -lngfem -lngcomp -lngsolve -lngla -o $@

clean:
	rm *.o libmyngsolve.so
