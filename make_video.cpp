#define FRAMERATE 20
#define STEPS_PER_FRAME 2
#define LENGTH 20

#define MACSIZE 150

#define IMG_H 600
#define IMG_W 600

#include <iostream>
#include <QPainter>

#include "mac_grid.hpp"
#include "colourUtil.hpp"

using std::cout;
using std::cerr;
using std::endl;

using std::string;

int main(int argc, char **argv)
{
    static int count = 0;

    if(argc != 2)
    {
        cerr << "Need to input directory name" << endl;
        return 1;
    }

    string prefix = argv[1];
    MacGrid mac(MACSIZE, 1.0f / (FRAMERATE * STEPS_PER_FRAME), true);

    for(int i = 0; i < LENGTH; ++i)
    {
        for(int j = 0; j < FRAMERATE; ++j)
        {
            QImage img(IMG_W, IMG_H, QImage::Format_RGB32);

            for(int k = 0; k < STEPS_PER_FRAME; ++k)
            {
                mac.step();
            }

            for(int i = 0; i < IMG_W; ++i)
            {
                for(int j = 0; j < IMG_H; ++j)
                {
                    int grid_width = ceil(IMG_W / (float)MACSIZE);
                    int grid_height = ceil(IMG_H / (float)MACSIZE);

                    uint colour;


                    int a = i / grid_width;
                    int b = j / grid_height;

                    if(!mac.isFuel(a, b) && mac.getTemp(a, b) > IGNITION_TEMP)
                    {
                        colour = getBlackbodyRGB( mac.getTemp(
                                    i / grid_width,
                                    j / grid_height) * 2.0); 
                        /* *
                           std::max(0.0f, mac.getDensity(
                           i / grid_width,
                           j / grid_height));*/
                    }
                    else if(!mac.isFuel(a, b))
                    {
                        colour = 255 * std::max(0.0f, mac.getDensity(
                                    i / grid_width,
                                    j / grid_height));
                        colour = colour + (colour << 8) + (colour << 16);
                    }
                    else
                    {
                        colour = 255;
                    }

                    img.setPixel(i, IMG_H - j - 1, colour);
                }
            }

            img.save(QString::fromStdString(prefix + "/frame" + std::to_string(count) + ".jpg"));
            ++count;
        }
    }
}
