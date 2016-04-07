#define FRAMERATE 10
#define STEPS_PER_FRAME 2
#define LENGTH 10

#define MACSIZE 100
#define SCALE 1.0f

#define IMG_H 600
#define IMG_W 600

#include <iostream>
#include <QPainter>

#include "mac_box.hpp"
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
    MacBox mac(MACSIZE, 1.0f / (FRAMERATE * STEPS_PER_FRAME), SCALE);

    for(int i = 0; i < LENGTH; ++i)
    {
        for(int j = 0; j < FRAMERATE; ++j)
        {
            for(int k = 0; k < STEPS_PER_FRAME; ++k)
            {
                mac.step();
            }
            
            QImage img(IMG_W, IMG_H, QImage::Format_RGB32);
            mac.computeImage(img);

            img.save(QString::fromStdString(prefix + "/frame" + std::to_string(count) + ".jpg"));
            cout << "Saving frame " << count << endl;
            ++count;
        }
    }
}
