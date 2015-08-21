
From PLUG mailing list:

> import -window root root.jpeg
> # will dump the screenshot to a jpeg file

Also, see "man import"

import -window root foo.jpeg

	Will capture the whole screen

import foo.jpeg

	Will wait for you do create a rectangle on the screen by holding down
	the mouse and dragging, and will then capture the rectangle when you
	release the mouse.

import -frame foo.jpeg

	Appears to do a good job of saving the window that you then click on.
	It might be necessary to add the "-border" parameter, but I haven't
	needed to do that yet.

