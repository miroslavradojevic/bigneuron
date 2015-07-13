// This macro imports a series of text images as a stack.

  dir = getDirectory("Choose directory");
  list = getFileList(dir);
print("there were " + list.length);
  run("Close All");
  setBatchMode(true);
  for (i=0; i<list.length; i++) {
     file = dir + list[i];
     run("Text Image... ", "open=&file");
  }
  run("Images to Stack", "use");
  setBatchMode(false);