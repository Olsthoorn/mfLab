classdef Colors
   properties
      R = 0;
      G = 0;
      B = 0;
   end
   methods
      function c = Colors(r, g, b)
         c.R = r; c.G = g; c.B = b;
      end
   end
   enumeration
      Blueish   (18/255,104/255,179/255)
      Reddish   (237/255,36/255,38/255)
      Greenish  (155/255,190/255,61/255)
      Purplish  (123/255,45/255,116/255)
      Yellowish (1,199/255,0)
      LightBlue (77/255,190/255,238/255)
   end
end