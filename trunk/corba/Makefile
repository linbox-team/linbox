LDFLAGS_SERVER = `orbit-config --libs server` \
		 `gnome-config --libs gnomeui gnorba`
LDFLAGS_CLIENT = `orbit-config --libs client` \
		 `gnome-config --libs gnomeui gnorba`
CFLAGS = `orbit-config --cflags client` `gnome-config --cflags gnomeui gnorba`

all: mymodule-client mymodule-server

mymodule-client: mymodule-client.o mymodule-common.o mymodule-stubs.o
	gcc $(LDFLAGS_SERVER) -o $@ $^

mymodule-server: mymodule-skelimpl.o mymodule-common.o mymodule-skels.o
	gcc $(LDFLAGS_CLIENT) -o $@ $^

mymodule-client.o: mymodule-client.c mymodule.h
	gcc $(CFLAGS) -c -o $@ $<

mymodule-common.o: mymodule-common.c mymodule.h
	gcc $(CFLAGS) -c -o $@ $<

mymodule-stubs.o: mymodule-stubs.c mymodule.h
	gcc $(CFLAGS) -c -o $@ $<

mymodule-skels.o: mymodule-skels.c mymodule.h
	gcc $(CFLAGS) -c -o $@ $<

mymodule-skelimpl.o: mymodule-skelimpl.c mymodule.h
	gcc $(CFLAGS) -c -o $@ $<

mymodule-skels.c mymodule-stubs.c mymodule-common.c mymodule.h: mymodule.idl
	orbit-idl mymodule.idl

clean:
	rm -f *.o mymodule-client mymodule-server
	rm -f mymodule-stubs.c mymodule-skels.c mymodule-common.c mymodule.h
